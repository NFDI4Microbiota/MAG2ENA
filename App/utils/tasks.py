import time
import streamlit as st
from enum import Enum
from pathlib import Path
from redis import Redis
from rq import Queue
from .db import TaskResultManager, TaskStatus
import requests
import os

redis_conn = Redis(host="localhost", port=6379)
q = Queue("ena", connection=redis_conn)

_APP_DIR = Path(__file__).resolve().parent.parent
manager = TaskResultManager(str(_APP_DIR / "task-results" / "task_results.db"))

_SENSITIVE_SUBMISSION_KEYS = {"ena_user", "ena_password"}
_SENSITIVE_JOB_KEYS = {"encryption_password"}

class JobStatus(Enum):
    PENDING = "pending"
    FINISHED = "finished"
    FAILED = "failed"
    INVALID = "invalid"


def scrub_job_credentials(job):
    """Remove ENA credentials from a persisted RQ job payload."""
    if job is None:
        return

    try:
        job_kwargs = dict(job.kwargs)
        for key in _SENSITIVE_JOB_KEYS:
            job_kwargs.pop(key, None)

        submission = job_kwargs.get("submission")
        if isinstance(submission, dict):
            sanitized_submission = {
                key: value
                for key, value in submission.items()
                if key not in _SENSITIVE_SUBMISSION_KEYS
            }
            job_kwargs["submission"] = sanitized_submission

        job.kwargs = job_kwargs
        job.save(include_result=False)
    except Exception:
        # Never log the job payload while handling credential cleanup.
        return

def send_job_email(email, job_id, status):
    """
    Send a job completion email using Gmail SMTP.
    """

    api_key = os.environ.get("API_KEY")

    # Email content
    if status == "started":
        subject = f"[MAG2ENA] Job submitted to queue"
        body = f"""Dear user,\n\nYour job was submitted to the queue.\nJob ID: {job_id}\n\nYou will receive an email when your job is finished."""
    elif status == "success":
        subject = f"[MAG2ENA] Job has finished"
        body = f"""Dear user,\n\nYour job has completed successfully.\nJob ID: {job_id}\n\nConsult the output in the Jobs module using your Job ID and, if encrypted, your results password."""
    elif status == "failed":
        subject = f"[MAG2ENA] Job has finished"
        body = f"""Dear user,\n\nYour job has failed.\n\nTry submitting again later."""

    response =  requests.post(
                "https://api.mailgun.net/v3/bioauto.inteligentehub.com.br/messages",
                auth=("api", api_key),
                data={"from": "MAG2ENA <submission@bioauto.inteligentehub.com.br>",
                    "to": email,
                    "subject": subject,
                    "text": body})
    
    if response.status_code == 200:
        # Print a message to console after successfully sending the email.
        print(f"Email sent to {email}.")
    else:
        print(f"Email failed to be sent to {email}.")

def _on_success(job, connection, result, *args, **kwargs):
    """
    RQ success callback signature: (job, connection, result)
    You can access the job kwargs with `job.kwargs` (a dict).
    """

    try:
        email = job.kwargs.get("email")
    except Exception:
        email = None

    scrub_job_credentials(job)
    manager.store_result(job.id, TaskStatus.SUCCESS)

    if email:
        send_job_email(email, job.id, "success")

def _on_failure(job, connection, *args, **kwargs):
    """
    RQ failure callback signature can be (job, *args...) — safer to ignore 'result' param name.
    """
    try:
        email = job.kwargs.get("email")
    except Exception:
        email = None

    scrub_job_credentials(job)
    manager.store_result(job.id, TaskStatus.FAILURE)

    if email:
        send_job_email(email, job.id, "failed")

def enqueue_task(fn, fn_kwargs=None):
    """
    Enqueue fn with fn_kwargs (dict). Always attach callbacks and
    always store a pending task entry in the DB.
    """
    if fn_kwargs is None:
        fn_kwargs = {}

    enqueue_kwargs = {
        "on_success": _on_success,
        "on_failure": _on_failure,
        # Avoid RQ's default description, which includes a repr of all arguments.
        "description": f"{fn.__module__}.{fn.__name__}",
        # Bound how long credentials can remain in Redis if no worker starts the job.
        "ttl": 86400,
        # Also bound retention if a job fails before its worker can scrub the payload.
        "failure_ttl": 86400,
    }

    # create the job
    job = q.enqueue(fn, kwargs=fn_kwargs, **enqueue_kwargs, job_timeout=7200)

    try:
        email = job.kwargs.get("email")
    except Exception:
        email = None

    if email:
        send_job_email(email, job.id, "started")
    
    id_ = job.id
    manager.store_pending_task(id_)

    return id_

def check_job_status(job_id):
    job = q.fetch_job(job_id)
    if job is None:
        return JobStatus.INVALID, None
    
    if job.is_finished:
        scrub_job_credentials(job)
        return JobStatus.FINISHED, job.result
    elif job.is_failed:
        # A worker import failure can occur before the normal failure callback
        # has a chance to scrub credentials. Terminal jobs no longer need them.
        scrub_job_credentials(job)
        return JobStatus.FAILED, job.exc_info
    else:
        return JobStatus.PENDING, None
