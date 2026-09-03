"""Password-based authenticated encryption for completed job artifacts."""

import base64
import io
import os
import secrets
import shutil
import tarfile
import tempfile
from pathlib import Path

from cryptography.fernet import Fernet, InvalidToken
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC


ARCHIVE_FILENAME = "job_archive.enc"
SALT_FILENAME = "job_salt.bin"
PBKDF2_ITERATIONS = 390_000


def derive_key_from_password(password: str, salt: bytes) -> bytes:
    """Derive a Fernet key without storing the user's password."""
    kdf = PBKDF2HMAC(
        algorithm=hashes.SHA256(),
        length=32,
        salt=salt,
        iterations=PBKDF2_ITERATIONS,
    )
    return base64.urlsafe_b64encode(kdf.derive(password.encode("utf-8")))


def is_encrypted_job(job_path: Path | str) -> bool:
    path = Path(job_path)
    return (path / ARCHIVE_FILENAME).is_file() and (path / SALT_FILENAME).is_file()


def _make_archive_bytes(job_path: Path) -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        for child in sorted(job_path.iterdir()):
            archive.add(child, arcname=child.name)
    return buffer.getvalue()


def _write_private_atomic(path: Path, content: bytes) -> None:
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=path.parent,
            prefix=f".{path.name}.",
            delete=False,
        ) as temp_file:
            temp_path = Path(temp_file.name)
            temp_path.chmod(0o600)
            temp_file.write(content)
            temp_file.flush()
            os.fsync(temp_file.fileno())
        os.replace(temp_path, path)
    finally:
        if temp_path:
            temp_path.unlink(missing_ok=True)


def _remove_path(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def encrypt_job_folder(
    job_path: Path | str,
    password: str,
    *,
    allowed_root: Path | str,
) -> None:
    """Encrypt all job artifacts and remove their plaintext originals.

    Encryption fails closed: if archive creation cannot finish, any plaintext
    artifacts in the job directory are removed before the error is propagated.
    """
    path = Path(job_path).resolve()
    root = Path(allowed_root).resolve()
    if not password:
        raise ValueError("An encryption password is required.")
    if not path.is_dir() or path.parent != root:
        raise ValueError("Invalid job directory.")

    archive_path = path / ARCHIVE_FILENAME
    salt_path = path / SALT_FILENAME

    try:
        archive_bytes = _make_archive_bytes(path)
        salt = secrets.token_bytes(16)
        encrypted = Fernet(derive_key_from_password(password, salt)).encrypt(
            archive_bytes
        )

        _write_private_atomic(archive_path, encrypted)
        _write_private_atomic(salt_path, salt)

        for child in list(path.iterdir()):
            if child not in {archive_path, salt_path}:
                _remove_path(child)
    except Exception:
        # Never leave unencrypted results behind for a password-protected job.
        for child in list(path.iterdir()):
            _remove_path(child)
        raise


def decrypt_job_archive(
    job_path: Path | str,
    password: str,
    target_path: Path | str,
) -> bool:
    """Decrypt a job archive into a private temporary directory."""
    path = Path(job_path)
    target = Path(target_path).resolve()

    try:
        salt = (path / SALT_FILENAME).read_bytes()
        encrypted = (path / ARCHIVE_FILENAME).read_bytes()
        archive_bytes = Fernet(
            derive_key_from_password(password, salt)
        ).decrypt(encrypted)

        target.mkdir(parents=True, exist_ok=True, mode=0o700)
        with tarfile.open(fileobj=io.BytesIO(archive_bytes), mode="r:gz") as archive:
            for member in archive.getmembers():
                member_path = (target / member.name).resolve()
                if target not in member_path.parents and member_path != target:
                    raise ValueError("Unsafe path in encrypted job archive.")
                if member.issym() or member.islnk():
                    raise ValueError("Links are not allowed in encrypted job archives.")
            archive.extractall(target, filter="data")
        return True
    except (InvalidToken, OSError, ValueError, tarfile.TarError):
        shutil.rmtree(target, ignore_errors=True)
        return False
