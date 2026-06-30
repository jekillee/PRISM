"""
PRISM headless batch CLI — invoked through the `prism` command.

Wired into run_prism.sh as the `nmode` subcommand, so on nkstar/ukstar you run:

    prism nmode --shot 40848 --tmin 2.0 --tmax 8.0
    prism nmode --shot 40848                                # full shot (omit --tmin/--tmax)
    prism nmode --shots 40848 40850 40852 --tmin 2.0 --tmax 8.0 --fmax 100
    prism nmode --shot-range 10000 10100                    # inclusive shot range

Each job computes an n-mode spectrum and stores it as an .npz file under the archive
root (same layout as the GUI's saved NPZ). When a shot's file already exists you are
asked whether to overwrite or skip it (with all / skip-all / quit for batches). One
failed shot never stops the batch.

This module never imports PySide6 — it is headless-only.
"""

import argparse
import math
import os
import sys

# Quiet MDSplus, skip the GUI's X11 setup (run_prism.sh still exports X11 vars,
# which are harmless here).
os.environ.setdefault("MDSPLUS_DEBUG", "0")

try:
    from config.app_config import PRISM_RESULTS_ROOT
except Exception:  # pragma: no cover - config should always import (Qt-free)
    PRISM_RESULTS_ROOT = os.path.expanduser("~/prism_results")


def _resolve_root():
    """Archive root: $PRISM_ARCHIVE_ROOT if set, else ~/prism_results."""
    return os.environ.get("PRISM_ARCHIVE_ROOT") or PRISM_RESULTS_ROOT


def _build_nmode_parser(subparsers):
    p = subparsers.add_parser(
        "nmode", help="Compute n-mode spectrum (one or many shots) and cache it"
    )
    p.add_argument("--shot", type=int, help="single shot number")
    p.add_argument("--shots", type=int, nargs="+", metavar="N",
                   help="multiple shot numbers")
    p.add_argument("--shot-range", type=int, nargs=2, metavar=("START", "END"),
                   help="inclusive shot range (e.g. --shot-range 10000 10100)")
    p.add_argument("--tmin", type=float, default=None,
                   help="window start [s] (default: full shot)")
    p.add_argument("--tmax", type=float, default=None,
                   help="window end [s] (default: full shot)")
    p.add_argument("--t-interval", type=float, default=0.01,
                   help="FFT window length [s] (default 0.01)")
    p.add_argument("--fmin", type=float, default=0.0, help="min frequency [kHz]")
    p.add_argument("--fmax", type=float, default=100.0, help="max frequency [kHz]")
    p.add_argument("--tol", type=float, default=0.8, help="residual tolerance")
    p.add_argument("--nmodes", type=int, default=5, help="number of modes")
    p.add_argument("--frac", type=float, default=1e-2, help="amplitude threshold fraction")
    p.add_argument("--msign", type=int, default=1, choices=[0, 1, -1, 2],
                   help="0=abs, 1=pos, -1=neg, 2=all (default 1)")
    p.add_argument("--integrate", action="store_true", help="integrate dB/dt -> B")
    p.add_argument("--no-detrend", dest="detrend", action="store_false",
                   help="disable per-window detrending")
    p.set_defaults(detrend=True)
    return p


def _shots_from_args(args):
    shots = []
    if args.shot is not None:
        shots.append(args.shot)
    if args.shots:
        shots.extend(args.shots)
    if args.shot_range:
        start, end = args.shot_range
        shots.extend(range(start, end + 1))   # END inclusive
    # de-duplicate, preserve order
    seen = set()
    return [s for s in shots if not (s in seen or seen.add(s))]


def _specs_from_args(args, shots):
    from batch import NModeJobSpec
    return [
        NModeJobSpec(
            shot=s, tmin=args.tmin, tmax=args.tmax, t_interval=args.t_interval,
            fmin=args.fmin, fmax=args.fmax, tol=args.tol, nmodes=args.nmodes,
            frac=args.frac, msign=args.msign, integrate=args.integrate,
            detrend=args.detrend,
        )
        for s in shots
    ]


def _fmt_shots(specs, cap=40):
    """Space-joined shot numbers, capped so a huge range doesn't flood the screen."""
    nums = [str(s.shot) for s in specs]
    if len(nums) <= cap:
        return " ".join(nums)
    return " ".join(nums[:cap]) + f" ... (+{len(nums) - cap} more)"


def _ask_overwrite(shot, path):
    """Prompt about an existing file. Returns 'yes' | 'no' | 'all' | 'skip-all' |
    'quit'. Non-interactive input (EOF) defaults to 'no' so a piped/cron run never
    clobbers existing files unattended."""
    prompt = (f"  shot {shot}: {path.name} exists. overwrite? "
              f"[y]es / [n]o / [a]ll / [s]kip-all / [q]uit: ")
    while True:
        try:
            ans = input(prompt).strip().lower()
        except EOFError:
            print()
            return "no"
        if ans in ("y", "yes"):
            return "yes"
        if ans in ("n", "no", ""):
            return "no"
        if ans in ("a", "all"):
            return "all"
        if ans in ("s", "skip", "skip-all"):
            return "skip-all"
        if ans in ("q", "quit"):
            return "quit"
        print("  please answer y / n / a / s / q")


def _run_nmode(args):
    shots = _shots_from_args(args)
    if not shots:
        print("error: provide --shot, --shots, or --shot-range", file=sys.stderr)
        return 2
    nonfinite = [n for n, v in (
        ("--tmin", args.tmin), ("--tmax", args.tmax), ("--t-interval", args.t_interval),
        ("--fmin", args.fmin), ("--fmax", args.fmax), ("--tol", args.tol),
        ("--frac", args.frac),
    ) if v is not None and not math.isfinite(v)]
    if nonfinite:
        print(f"error: non-finite value for {', '.join(nonfinite)}", file=sys.stderr)
        return 2

    specs = _specs_from_args(args, shots)
    root = _resolve_root()
    from batch.archive import archive_path

    # Categorize first so the user sees which shots already exist before deciding.
    new_specs, existing = [], []
    for spec in specs:
        (existing if archive_path(root, spec).exists() else new_specs).append(spec)

    if new_specs:
        print(f"{len(new_specs)} new shot(s) to compute: {_fmt_shots(new_specs)}")
    if existing:
        print(f"{len(existing)} shot(s) already exist: {_fmt_shots(existing)}")

    # New shots compute outright; decide each existing shot (overwrite/skip).
    to_run, skipped = list(new_specs), []
    bulk = None   # None -> ask each; "all" -> overwrite rest; "skip-all" -> skip rest
    for spec in existing:
        if bulk == "all":
            to_run.append(spec)
            continue
        if bulk == "skip-all":
            skipped.append(spec)
            continue
        ans = _ask_overwrite(spec.shot, archive_path(root, spec))
        if ans == "quit":
            print("aborted; nothing computed.", file=sys.stderr)
            return 1
        if ans == "all":
            bulk = "all"
            to_run.append(spec)
        elif ans == "skip-all":
            bulk = "skip-all"
            skipped.append(spec)
        elif ans == "yes":
            to_run.append(spec)
        else:  # "no"
            skipped.append(spec)

    for spec in skipped:
        print(f"[skipped]  shot {spec.shot}")

    if not to_run:
        print(f"\nNothing computed: {len(skipped)} skipped  ->  {root}")
        return 0

    from batch import run_many

    def _progress(outcome):
        if outcome.status == "failed":
            print(f"[failed]   shot {outcome.spec.shot}  ({outcome.seconds:.1f}s)  "
                  f"{outcome.error}")
        else:
            print(f"[computed] shot {outcome.spec.shot}  ({outcome.seconds:.1f}s)")

    outcomes = run_many(
        to_run, archive_root=root, overwrite=True, progress=_progress,
    )
    n_ok = sum(o.status in ("computed", "cached") for o in outcomes)
    n_fail = sum(o.status == "failed" for o in outcomes)
    print(f"\nDone: {n_ok} computed, {len(skipped)} skipped, {n_fail} failed  ->  {root}")
    return 1 if n_fail and not n_ok else 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="prism", description="PRISM headless batch CLI")
    subparsers = parser.add_subparsers(dest="cmd", required=True)
    _build_nmode_parser(subparsers)

    args = parser.parse_args(argv)
    if args.cmd == "nmode":
        return _run_nmode(args)
    parser.error(f"unknown command: {args.cmd}")


if __name__ == "__main__":
    sys.exit(main())
