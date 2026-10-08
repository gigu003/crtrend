# Native multi-cohort permutation reference

Imported on 2026-10-07 from actual V6.1.0 text exports under
`/Users/qc0824/Documents/IMS/Joinpoint/Selection-audit-results/Perm-groups`:

- Original/groups-original.Export.Model.Selection.txt
- relabelled/groups-relabelled.Export.Model.Selection.txt
- single-0/single-0.Export.Model.Selection.txt
- single-1/single-1.Export.Model.Selection.txt
- single-2/single-2.Export.Model.Selection.txt

Each prefix also has Run.Summary, Final.Selected.Model and Data text checked
by data-raw/joinpoint-permutation-groups-reference.R. The fixture contains
all 36 selection rows and associated final model counts, not R-generated P values.
Original cohort=sex; relabelled sex=2-cohort; single-N sex=N.
Settings: provided age-adjusted rates and SE, log transform, independent errors,
grid 0–4 joinpoints, minimum end/between observations 2, no inserted grid points,
fixed 4499 permutations, seed7160, alpha0.05. Export count4500 includes original.
Only this independent-error, fixed-seed configuration is externally benchmarked;
tSPRT and correlated errors are not covered.
