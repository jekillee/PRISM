"""
Run gate for the PRISM batch SDK.

`run` executes a single job cache-first: it returns a cached result if one exists
at the spec's archive path, otherwise computes, stores, and returns it.
`run_many` drives a sequential batch where one failed shot never aborts the rest.

Jobs run sequentially on purpose: load_mirnov_data already opens one MDS+
connection per Mirnov channel, so concurrent shots could flood the MDS+ server.
"""

import os
import time as _time
from dataclasses import dataclass

from batch.archive import get_or_compute
from batch.compute import compute_nmode
from batch.results import NModeResult

try:
    from config.app_config import PRISM_TMP_ROOT, ensure_shared_dir
except Exception:  # pragma: no cover - config should always import (Qt-free)
    PRISM_TMP_ROOT = "/tmp/prism"
    def ensure_shared_dir(path, mode=0o777):
        os.makedirs(path, exist_ok=True)
        return path

# subsystem -> (compute_fn(spec, *, mds_server), result_cls)
_DISPATCH = {
    "nmode": (compute_nmode, NModeResult),
}


@dataclass
class JobOutcome:
    spec: object
    status: str            # "cached" | "computed" | "failed"
    seconds: float = 0.0
    error: str = ""
    result: object = None


def _resolve_archive_root(archive_root):
    # Precedence: explicit arg > $PRISM_ARCHIVE_ROOT > the shared /tmp/prism default.
    root = archive_root or os.environ.get("PRISM_ARCHIVE_ROOT") or PRISM_TMP_ROOT
    # When using the shared default tree, make the root + nmode/ subdir world-writable
    # so multiple users on nkstar/ukstar share one /tmp/prism (nmode files land in
    # <root>/nmode/<shot>/). A user-supplied root is left untouched.
    if os.path.abspath(root) == os.path.abspath(PRISM_TMP_ROOT):
        ensure_shared_dir(PRISM_TMP_ROOT)
        ensure_shared_dir(os.path.join(root, "nmode"))
    return root


def run(spec, *, archive_root=None, mds_server=None, overwrite=False):
    """Compute or load one job. Returns (result, status) with status
    'cached' (existing file reused) or 'computed'. overwrite=True recomputes."""
    root = _resolve_archive_root(archive_root)
    try:
        compute_fn, result_cls = _DISPATCH[spec.SUBSYSTEM]
    except KeyError:
        raise ValueError(f"No compute registered for subsystem '{spec.SUBSYSTEM}'")

    def _fn(s):
        return compute_fn(s, mds_server=mds_server)

    return get_or_compute(root, spec, _fn, result_cls, overwrite=overwrite)


def run_many(specs, *, archive_root=None, mds_server=None, overwrite=False,
             on_error="continue", progress=None):
    """Run many jobs sequentially. Returns a list of JobOutcome.

    With on_error='continue' (default), a failing job is recorded as
    status='failed' and the batch proceeds; otherwise the exception propagates.
    `progress` is an optional callback invoked with each JobOutcome.
    """
    root = _resolve_archive_root(archive_root)
    outcomes = []
    for spec in specs:
        t0 = _time.time()
        try:
            result, status = run(spec, archive_root=root,
                                  mds_server=mds_server, overwrite=overwrite)
            outcome = JobOutcome(spec, status, _time.time() - t0, result=result)
        except Exception as exc:
            outcome = JobOutcome(spec, "failed", _time.time() - t0, error=repr(exc))
            if on_error != "continue":
                raise
        outcomes.append(outcome)
        if progress is not None:
            progress(outcome)
    return outcomes
