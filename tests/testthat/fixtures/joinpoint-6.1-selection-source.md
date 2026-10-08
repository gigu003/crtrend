# Source of joinpoint-6.1-selection.csv

User-supplied Joinpoint 6.1.0 exports, received 2026-10-07:
`/Users/qc0824/Documents/IMS/Joinpoint/Selection-audit-results/`

- MBIC/input.Export.Model.Selection.txt: all five MBIC scores, seven decimals.
- DDS/input.Export.Selection.Results.txt: Delta BIC and Delta BIC3 both
  15.5408412, both selected k=1. No Delta reference for other k; stored as NA.
- DDS/input.Export.Run.Summary.txt and MBIC/input.Export.Run.Summary.txt:
  input.txt, years 2000-2021, provided rate/SE, log response, independent errors,
  grid search, min/max k=0/4, both min-observation settings=2, no extra grid
  points, no PPD. Input is jp_data sex=0. MBIC used EQ, DDS Parametric; this
  does not change the candidate selection scores being benchmarked here.

MBIC native score differs from the published expression by -2*log(n)/n,
matching the static native arithmetic at 0x4ca804-0x4ca817 (2*k log(n)/n).
The regression fixture checks production calculations, not adjusted reference
values. It establishes this case only; DDS threshold branches are not exercised
because BIC and BIC3 choose the same k.
# Permutation reference added 2026-10-07

`joinpoint-6.1-permutation.csv` is transcribed from
`/Users/qc0824/Documents/IMS/Joinpoint/Selection-audit-results/Perm/export.Export.Model.Selection.txt`.
The accompanying `export.Export.Run.Summary.txt` records V6.1.0,
JPO creation 2026-10-07 08:42:21, seed 7160, 4499 random permutations,
alpha 0.05, log rates, supplied SE, independent errors, grid 0–4 joinpoints,
minimum end/between observations 2, and no inserted grid points.
Input is the sex=0 subset of jp_data. Exported count 4500 includes original data.
Inference was EQ/5001, separate from these selection p-values.
Native stream and forward cumulative shuffle were traced at 0x4679d0,
0x4c7b40, 0x481e40, and 0x471980 in the previously hashed 6.1 binary.
No native compatibility claim is made for tSPRT or non-independent covariance
structures. Multi-group scheduling was subsequently benchmarked separately;
see joinpoint-6.1-permutation-groups-source.md.
