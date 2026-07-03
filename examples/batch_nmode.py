#!/usr/bin/env python3
"""Example: compute an n-mode spectrum with the PRISM batch API and inspect it.

Run on nkstar/ukstar (inside the KSTAR network) with the `prism python` helper, which
supplies PRISM's bundled Python + import path (no env setup needed):

    prism python examples/batch_nmode.py       # from the PRISM install directory

Results cache under ~/prism_results/nmode/ (override with $PRISM_ARCHIVE_ROOT).
"""
from batch import NModeJobSpec, run

# WHAT to compute: shot + time window + frequency band. Omit tmin/tmax for full shot.
spec = NModeJobSpec(shot=40848, tmin=2.0, tmax=8.0, fmin=0, fmax=100)

result, status = run(spec)
print(f"status = {status}\n")

# Printing the result summarizes every field and its shape.
print(result)

# amplitude_max[i] / amplitude_sum[i] correspond to n_modes_list[i].
nlist = result.meta["n_modes_list"]
print(f"\nmodes: {nlist}")
for i, n in enumerate(nlist):
    print(f"  n={n:+d}: peak Max amplitude = {result.amplitude_max[i].max():.3e} Gauss")
