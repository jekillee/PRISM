"""
PRISM headless batch SDK.

Run expensive diagnostics computations without the GUI and cache each result as an
.npz file (one per shot, same layout as the GUI's saved NPZ) so repeated requests
never recompute. Currently supports n-mode (toroidal mode-number) spectra; more
operations register into the same dispatch.

Quick start (run on nkstar/ukstar, inside the KSTAR network):

    from batch import NModeJobSpec, run

    spec = NModeJobSpec(shot=40848, tmin=2.0, tmax=8.0, fmin=0, fmax=100)
    result, status = run(spec, archive_root="/home/users/jklee/PRISM_archive")
    # status == "computed" the first time, "cached" afterwards
    # result is an NModeResult (time, frequency, mode_spectrum, amplitude, meta)

Sequential batch over many shots:

    from batch import NModeJobSpec, run_many

    specs = [NModeJobSpec(shot=s, tmin=2.0, tmax=8.0) for s in (40848, 40850, 40852)]
    outcomes = run_many(specs, archive_root="/home/users/jklee/PRISM_archive",
                        progress=lambda o: print(o.status, o.spec.shot))

The archive_root may also be supplied via the PRISM_ARCHIVE_ROOT environment
variable. Today it is a plain directory on each server; once the shared NAS is
mounted, point archive_root at the NAS mount (e.g. /mnt/prism) — no code change.
"""

from batch.jobspec import JobSpec, NModeJobSpec
from batch.results import NModeResult

__all__ = [
    "JobSpec",
    "NModeJobSpec",
    "NModeResult",
    "run",
    "run_many",
    "JobOutcome",
]


def __getattr__(name):
    # Lazy-load the run gate (PEP 562) so importing the package — `from batch import
    # NModeJobSpec`, spec building, --dry-run, cache-path listing — never pulls in the
    # compute stack (batch.runner -> batch.compute -> core.nmode -> scipy / MDSplus).
    # Only an actual run() / run_many() reference triggers those heavy imports.
    if name in ("run", "run_many", "JobOutcome"):
        from batch.runner import run, run_many, JobOutcome
        return {"run": run, "run_many": run_many, "JobOutcome": JobOutcome}[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
