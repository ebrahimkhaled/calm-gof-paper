# Reproduction archive

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22467285.svg)](https://doi.org/10.5281/zenodo.22467285)


**Paper:** A closed-form reference distribution for goodness-of-fit testing under penalised
logistic regression in the proportional regime.

Everything in the paper and the supplement is produced by the scripts here. There is no hidden
step: each table and figure names the file it comes from, and every simulation writes both a
per-replication file and a summary.

## What is here

| Directory | Contents |
|---|---|
| `R/` | `closedform_engine.R`, the reference implementation. `calm.gof()` is the shipped test; the file also carries the alternative references, groupings and inflations that the paper compares, because the case for the shipped design is made by comparison. |
| `sims/` | One script per experiment. The name is the label used in the paper's text and in `AI_BRAIN.md`. |
| `results/` | Every `*_perrep.csv` (one row per replication) and `*_summary.csv` (the aggregates quoted in the paper), plus the run logs. |
| `paper/` | `paper.tex`, `supplement.tex`, the figure scripts, and the audit script that checks the manuscript against the result files. |

## Reproducing

R 4.4 or later. Packages: `CompQuadForm` (exact tail), `ebrahim.gof` (the prepivoting comparator
`shrink.gof()` and the released `calm.gof()`), `GRPtests` (the rival test), `TH.data` (the glaucoma
data), `logistf` (checking the Firth fit), `parallel`.

```
Rscript sims/E4c3_final_size_glaucoma.R 500 8      # Table 1 and the glaucoma p-values
Rscript sims/E6_nongaussian.R           400 8      # Table 2
Rscript sims/E10_heldout_validation.R   500 300 8  # Table 3, and the sparse comparison
Rscript sims/E12_referee_items.R        500 8      # Table 4, the scale recovery, the tau study
Rscript sims/E4c_power_survivors.R      400 299 8  # Table 5 and Figure 4
Rscript sims/E4d_rivals.R               300 8      # Tables 5 and 6, the rival tests
Rscript sims/E4b2_calm_final_on_E4b.R   500 6      # Table 7
Rscript sims/E7_penalty_generality.R    300 8      # Table S1
Rscript paper/R/figures.R                          # Figures 1, 2, 3, 5
Rscript paper/R/fig5_power.R                       # Figure 4 and the power tables
python  paper/R/audit_seriesB.py                   # rebuilds the PDF and checks it
```

Every script fixes its own seed at the top; the seed is the same across scripts that must see the
same datasets, so paired comparisons are genuinely paired. Times are for eight workers on a
desktop: most scripts are minutes, the two that call the bootstrap are one to five hours.

`sessionInfo()` at the time the results were produced is in `results/sessionInfo.txt`.
