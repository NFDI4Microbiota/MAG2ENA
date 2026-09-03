import streamlit as st
import json
import time
import shutil
import tempfile
from pathlib import Path
from datetime import datetime
from utils.tasks import manager, check_job_status, JobStatus
from utils.db import TaskStatus
from utils.encryption import decrypt_job_archive, is_encrypted_job

_HERE = Path(__file__).parent
_APP_DIR = _HERE.parent
SUBMISSIONS_DIR = _APP_DIR / "jobs"
EXAMPLE_JOBS_DIR = _APP_DIR / "example_jobs"
EXAMPLE_JOB_ID = "a66917de-84d1-4751-b17c-a3c6657f6051"
_UNLOCKED_JOB_DIR_KEY = "_unlocked_job_dir"
_UNLOCKED_JOB_ID_KEY = "_unlocked_job_id"


def clear_decrypted_job():
    """Remove plaintext artifacts extracted for the current browser session."""
    decrypted_path = st.session_state.pop(_UNLOCKED_JOB_DIR_KEY, None)
    st.session_state.pop(_UNLOCKED_JOB_ID_KEY, None)
    if not decrypted_path:
        return

    path = Path(decrypted_path).resolve()
    temp_root = Path(tempfile.gettempdir()).resolve()
    if temp_root in path.parents and path.name.startswith("mag2ena_job_"):
        shutil.rmtree(path, ignore_errors=True)


def _unlocked_job_dir(job_id):
    if st.session_state.get(_UNLOCKED_JOB_ID_KEY) != job_id:
        return None
    path_value = st.session_state.get(_UNLOCKED_JOB_DIR_KEY)
    if not path_value:
        return None
    path = Path(path_value)
    return path if path.is_dir() else None


def _unlock_job(job_id, password):
    clear_decrypted_job()
    target = Path(tempfile.mkdtemp(prefix="mag2ena_job_"))
    target.chmod(0o700)
    if not decrypt_job_archive(_job_dir(job_id), password, target):
        shutil.rmtree(target, ignore_errors=True)
        return None

    st.session_state[_UNLOCKED_JOB_ID_KEY] = job_id
    st.session_state[_UNLOCKED_JOB_DIR_KEY] = str(target)
    return target


def _job_dir(job_id):
    if job_id == EXAMPLE_JOB_ID:
        return EXAMPLE_JOBS_DIR / job_id
    return SUBMISSIONS_DIR / job_id


def _format_ts(ts_str):
    if not ts_str:
        return "—"
    try:
        dt = datetime.fromisoformat(ts_str)
        return dt.strftime("%Y-%m-%d %H:%M UTC")
    except Exception:
        return ts_str


def _compute_duration(start_str, end_str):
    if not start_str or not end_str:
        return "—"
    try:
        start = datetime.fromisoformat(start_str)
        end = datetime.fromisoformat(end_str)
        total = int((end - start).total_seconds())
        h = total // 3600
        m = (total % 3600) // 60
        s = total % 60
        return f"{h:02d}:{m:02d}:{s:02d}"
    except Exception:
        return "—"


def _failure_message(content):
    """Return an actionable failure message without exposing a Python traceback."""
    text = str(content or "").strip()
    if not text:
        return "No additional error details were recorded."
    if "Traceback (most recent call last)" not in text:
        return text
    final_line = next(
        (line.strip() for line in reversed(text.splitlines()) if line.strip()), ""
    )
    if ": " in final_line:
        final_line = final_line.split(": ", 1)[1]
    return final_line or "The background worker encountered an unexpected error."


def _show_completed_job(job_id, db_info, output_dir=None):
    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric("Status", "Success")

    with col2:
        st.metric("Completed", _format_ts(db_info["end_time"] if db_info else None))

    with col3:

        duration = _compute_duration(
            db_info["start_time"] if db_info else None,
            db_info["end_time"] if db_info else None,
        )
        st.metric("Duration", duration)

    log_dir = Path(output_dir) if output_dir else _job_dir(job_id)
    result_file = log_dir / "result.json"

    if not result_file.exists():
        st.warning(
            "Result data not available. "
            "The submission output may have been cleaned up."
        )
        return

    with open(result_file) as f:
        result = json.load(f)

    job_type = result.get("job_type", "submission")
    col1, col2 = st.columns(2)

    with col1:
        st.metric(
            "Samples updated" if job_type == "metadata_update" else "Samples submitted",
            result.get("modified", result.get("submitted", 0)),
        )

    with col2:
        if job_type == "metadata_update":
            st.metric("Metadata cells changed", result.get("changed_cells", 0))
        else:
            st.metric("Errors", result.get("errors", 0))

    if not log_dir.exists():
        st.warning("Log directory not found.")
        return

    col_xml, col_zip = st.columns(2)

    submit_xml = log_dir / "submit.xml"
    if submit_xml.exists():
        with col_xml:
            st.info(
                "The XML document sent to ENA's submission API. "
                "Contains all sample metadata, study information, and checklist attributes.",
                icon="ℹ️",
            )
            st.download_button(
                label="Download submit.xml",
                data=submit_xml.read_bytes(),
                file_name="submit.xml",
                mime="text/xml",
                width="stretch",
                key=f"dl_{job_id}_submit.xml",
            )

    manifests_zip = log_dir / "manifests.zip"
    if manifests_zip.exists():
        with col_zip:
            st.info(
                "One Webin-CLI manifest file per MAG. "
                "Each manifest lists the study, sample accession, coverage, "
                "assembly software, sequencing platform, and FASTA path used for the assembly submission.",
                icon="ℹ️",
            )
            st.download_button(
                label="Download manifests (ZIP)",
                data=manifests_zip.read_bytes(),
                file_name="manifests.zip",
                mime="application/zip",
                width="stretch",
                key=f"dl_{job_id}_manifests.zip",
            )

    provenance_tsv = log_dir / "provenance.tsv"
    provenance_json = log_dir / "provenance.json"
    if provenance_tsv.exists() or provenance_json.exists():
        st.subheader("Provenance")
        st.info(
            "Unified history of metadata changes, resolver decisions, validation, "
            "file transformations, generated artifacts, and ENA submission outcomes. "
            "Use TSV for quick review or JSON for structured processing.",
            icon="ℹ️",
        )
        col_tsv, col_json = st.columns(2)
        if provenance_tsv.exists():
            with col_tsv:
                st.download_button(
                    label="Download provenance.tsv",
                    data=provenance_tsv.read_bytes(),
                    file_name="provenance.tsv",
                    mime="text/tab-separated-values",
                    width="stretch",
                    key=f"dl_{job_id}_provenance.tsv",
                )
        if provenance_json.exists():
            with col_json:
                st.download_button(
                    label="Download provenance.json",
                    data=provenance_json.read_bytes(),
                    file_name="provenance.json",
                    mime="application/json",
                    width="stretch",
                    key=f"dl_{job_id}_provenance.json",
                )

    _LOG_DESCRIPTIONS = {
        "webin_log.xml": (
            "ENA's response to the metadata-only MODIFY request, including validation messages."
            if job_type == "metadata_update"
            else "ENA's response to the metadata XML submission. Contains the accession "
            "numbers assigned to each sample and any processing messages."
        ),
        "success.txt": (
            "Webin-CLI output for assemblies that were submitted successfully. "
            "Each entry shows the sample alias and the tool's confirmation message."
        ),
        "error.txt": (
            "Error raised while retrieving, validating, or updating existing sample metadata."
            if job_type == "metadata_update"
            else "Webin-CLI output for assemblies that failed to submit. Check this file "
            "for validation errors or rejected fields."
        ),
        "before.xml": (
            "Authoritative ENA sample XML retrieved before the metadata update."
        ),
        "modify.xml": (
            "Metadata-only MODIFY document sent to ENA. It contains no assembly files or manifests."
        ),
        "rollback.xml": (
            "MODIFY document containing the original values. Review it carefully before any manual rollback."
        ),
        "metadata_diff.tsv": (
            "Field-level comparison of the original ENA values and curator-approved replacements."
        ),
        "validation_receipt.xml": (
            "ENA receipt from the non-mutating MODIFY + VALIDATE dry run."
        ),
    }

    log_files = [
        p for p in sorted(log_dir.glob("*"))
        if p.name not in {
            "result.json",
            "submit.xml",
            "manifests.zip",
            "job_info.json",
            "provenance.tsv",
            "provenance.json",
        }
    ]

    if log_files:
        st.subheader("Submission logs")

        for logfile in log_files:
            suffix = logfile.suffix.lower()

            try:
                content = logfile.read_text(encoding="utf-8")
            except Exception:
                content = None

            with st.expander(logfile.name, expanded=True):
                description = _LOG_DESCRIPTIONS.get(logfile.name)
                if description:
                    st.info(description, icon="ℹ️")

                if content:
                    lang = {".xml": "xml", ".json": "json"}.get(suffix, "text")
                    st.code(content, language=lang)
                else:
                    st.warning("Could not read file content.")

                with open(logfile, "rb") as f:
                    st.download_button(
                        label=f"Download {logfile.name}",
                        data=f.read(),
                        file_name=logfile.name,
                        width="stretch",
                        key=f"dl_{job_id}_{logfile.name}",
                    )


def _load_example_info():
    info_file = EXAMPLE_JOBS_DIR / EXAMPLE_JOB_ID / "job_info.json"
    with open(info_file) as f:
        return json.load(f)


def _show_job(job_id, password=None, unlock_requested=False):
    rq_details = None
    if job_id == EXAMPLE_JOB_ID:
        # Self-contained: doesn't rely on Redis or task_results.db.
        db_info = _load_example_info()
        status = JobStatus.FINISHED if db_info["status"] == "success" else JobStatus.FAILED
    else:
        status, rq_details = check_job_status(job_id)
        db_info = manager.get_result(job_id)

        # Neither Redis nor DB knows about this job
        if status == JobStatus.INVALID and db_info is None:
            st.error(f"Job `{job_id}` not found.")
            return

        # Redis is authoritative while an RQ job still exists. Reconcile cases
        # where a worker failed before its success/failure callback imported.
        if status == JobStatus.FAILED and (
            not db_info or db_info.get("status") not in {"success", "failure"}
        ):
            manager.store_result(job_id, TaskStatus.FAILURE)
            db_info = manager.get_result(job_id)
        elif status == JobStatus.FINISHED and (
            not db_info or db_info.get("status") not in {"success", "failure"}
        ):
            manager.store_result(job_id, TaskStatus.SUCCESS)
            db_info = manager.get_result(job_id)

    # Prefer DB status for completed jobs (Redis results expire)
    db_status = db_info["status"] if db_info else None
    terminal = (
        db_status in {"success", "failure"}
        or status in {JobStatus.FINISHED, JobStatus.FAILED}
    )
    stored_job_dir = _job_dir(job_id)
    output_dir = None

    if terminal and is_encrypted_job(stored_job_dir):
        output_dir = _unlocked_job_dir(job_id)

        if output_dir is None and unlock_requested:
            if not password:
                st.error("Enter the results password to unlock this encrypted job.")
                return
            output_dir = _unlock_job(job_id, password)
            if output_dir is None:
                st.error("Incorrect password or damaged encrypted job archive.")
                return

        if output_dir is None:
            st.info(
                "🔒 This job is encrypted. Enter its results password above and "
                "select **Open job** to view or download the results."
            )
            return

    if db_status == "success" or status == JobStatus.FINISHED:
        if output_dir:
            _show_completed_job(job_id, db_info, output_dir=output_dir)
        else:
            if st.session_state.get(_UNLOCKED_JOB_ID_KEY) != job_id:
                clear_decrypted_job()
            _show_completed_job(job_id, db_info)

    elif db_status == "failure" or status == JobStatus.FAILED:
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Status", "Failed")
        with col2:
            st.metric("Started", _format_ts(db_info["start_time"] if db_info else None))
        with col3:
            st.metric("Ended", _format_ts(db_info["end_time"] if db_info else None))
        st.error(
            "The submission job failed. "
            "Please check your credentials and metadata and try again."
        )

        error_file = (Path(output_dir) if output_dir else stored_job_dir) / "error.txt"
        if error_file.exists():
            try:
                content = error_file.read_text(encoding="utf-8")
            except Exception:
                content = None
            with st.expander("error.txt", expanded=True):
                safe_content = _failure_message(content)
                st.warning(safe_content)
                st.download_button(
                    label="Download error.txt",
                    data=safe_content.encode("utf-8"),
                    file_name="error.txt",
                    width="stretch",
                    key=f"dl_{job_id}_error.txt",
                )
        # Raw RQ/Python tracebacks are deliberately never rendered to users.

    else:
        pos = manager.get_job_position(job_id)
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Status", "Queued / Running")
        with col2:
            st.metric(
                "Queue position",
                str(pos) if pos is not None else "—",
            )
        st.info(
            "Your submission is being processed. "
            "This page will refresh automatically every 5 seconds."
        )
        time.sleep(5)
        st.rerun()

def runUI():

    if "job_input" not in st.session_state and "last_job_id" in st.session_state:
        st.session_state["job_input"] = st.session_state["last_job_id"]

    # st.subheader("Queue overview")

    # col1, col2 = st.columns(2)

    # with col1:
    #     st.markdown("**Pending & running**")
    #     pending_df = manager.get_pending_jobs()
    #     if pending_df.empty:
    #         st.info("No jobs currently in queue.")
    #     else:
    #         st.dataframe(
    #             pending_df,
    #             hide_index=True,
    #             width="stretch",
    #         )

    # with col2:
    #     st.markdown("**Recently completed**")
    #     completed_df = manager.get_recent_completed_jobs()
    #     if completed_df.empty:
    #         st.info("No completed jobs yet.")
    #     else:
    #         st.dataframe(
    #             completed_df,
    #             hide_index=True,
    #             width="stretch",
    #         )

    # st.divider()

    # st.subheader("Look up a submission")

    def _set_example():
        st.session_state["job_input"] = EXAMPLE_JOB_ID

    with st.container(border=True):
        col1, col2 = st.columns([9, 1])

        with col2:
            st.button(
                "Example",
                width="stretch",
                on_click=_set_example,
            )

        with st.form("jobs_submit", border=False, clear_on_submit=True):
            col_job, col_password = st.columns(2)
            with col_job:
                job_id = st.text_input("Enter Job ID", key="job_input")
            with col_password:
                password = st.text_input(
                    "Results password (if encrypted)",
                    type="password",
                    help=(
                        "Enter the password created during submission. Unencrypted jobs "
                        "do not require a password."
                    ),
                )
            submitted = st.form_submit_button(
                "Open job",
                width="stretch",
                type="primary",
            )

    if submitted and job_id:
        selected_job_id = job_id.strip()
        if st.session_state.get(_UNLOCKED_JOB_ID_KEY) != selected_job_id:
            clear_decrypted_job()
        st.session_state["job_id"] = selected_job_id

    if "job_id" in st.session_state:
        _show_job(
            st.session_state["job_id"],
            password=password if submitted else None,
            unlock_requested=submitted,
        )
