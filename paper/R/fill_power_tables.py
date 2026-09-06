# Fill Table tab:power (valid tests, raw rates) and Table tab:invalid (size + size-adjusted power) from results/power_master_summary.csv
# (written by paper/R/fig5_power.R).  Replaces the two "%% PENDING E4d rows" placeholders and removes the other PENDING E4d markers.
import io, csv, re, os, sys
P = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "paper.tex")
RES = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results", "power_master_summary.csv")
if len(sys.argv) > 1: P = sys.argv[1]          # dry-run on a copy: python fill_power_tables.py <copy-of-paper.tex> [<summary.csv>]
if len(sys.argv) > 2: RES = sys.argv[2]
rows = list(csv.DictReader(io.open(RES, encoding="utf-8")))
cells = ["fixed p (κ = 0.01)", "κ = 0.10", "κ = 0.25", "κ = 0.40"]
cell_lab = {"fixed p (κ = 0.01)": "A ($\\kappa = 0.01$)", "κ = 0.10": "B, $\\kappa = 0.10$", "κ = 0.25": "B, $\\kappa = 0.25$", "κ = 0.40": "B, $\\kappa = 0.40$"}
dep_lab = {"null": "---", "quadratic": "quadratic in the index", "offindex": "$x_1x_2$ off the index"}
order = {"null": 0, "quadratic": 1, "offindex": 2}
def num(r, k):
    v = r.get(k, "")
    try: return float(v)
    except: return None
def f(x): return "---" if x is None else "%.3f" % x
# validity per (cell, test): size in that cell <= .075 ; MLE column also requires the MLE to exist in >= 80% of datasets
size = {(r["cell"], t): num(r, t) for r in rows if r["dep"] == "null" for t in ["calm_hl", "calm_edge", "pp_hl", "pp_edge", "grp", "mle_hl", "mle_edge", "cvhl", "calslope", "spieg", "uss", "unc_hl"]}
mle_ok = {r["cell"]: num(r, "mle_ok") for r in rows if r["dep"] == "null"}
def cellv(c, t):
    s = size.get((c, t)); ok = s is not None and s <= .075
    if t.startswith("mle"):   # the two MLE columns stand or fall together: the fit must exist and BOTH sizes must hold
        s2 = [size.get((c, "mle_hl")), size.get((c, "mle_edge"))]
        ok = all(v is not None and v <= .075 for v in s2) and (mle_ok.get(c) or 0) >= .8
    return ok
valid_cols = ["calm_hl", "calm_edge", "pp_hl", "pp_edge", "grp", "mle_hl", "mle_edge"]
lines = []
for c in cells:
    rs = sorted([r for r in rows if r["cell"] == c], key=lambda r: (order[r["dep"]], float(r["rel"])))
    for i, r in enumerate(rs):
        lab = cell_lab[c] if i == 0 else ""
        rel = "0" if r["dep"] == "null" else ("%g" % float(r["rel"]))
        vals = [f(num(r, t)) if (cellv(c, t) or r["dep"] == "null") else "---" for t in valid_cols]
        lines.append("%s & %s & %s & %s \\\\" % (lab, dep_lab[r["dep"]], rel, " & ".join(vals)))
    if c != cells[-1]: lines.append("\\addlinespace")
power_tab = ("\\begin{tabular}{lllccccccc}\n\\toprule\n& & & \\multicolumn{2}{c}{CALM} & \\multicolumn{2}{c}{prepivot} & & \\multicolumn{2}{c}{MLE fit} \\\\\n"
             "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){9-10}\nDesign & departure & size & SC.HL & SC.EDGE & SC.HL & SC.EDGE & GRP & HL & EDGE \\\\\n\\midrule\n"
             + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}")
# invalid table: size (null) and size-adjusted power at quadratic rel 1, per cell
inv_cols = [("unc_hl", "uncorrected Hosmer--Lemeshow"), ("calslope", "calibration intercept and slope"), ("spieg", "Spiegelhalter's $z$"), ("uss", "unweighted sum of squares"), ("cvhl", "Hosmer--Lemeshow on cross-validated predictions"), ("mle_hl", "Hosmer--Lemeshow on the MLE fit")]
head = " & ".join(["\\multicolumn{2}{c}{%s}" % ("A" if c.startswith("fixed") else "$\\kappa = %s$" % c.split("= ")[1]) for c in cells])
il = []
for t, lab in inv_cols:
    vals = []
    for c in cells:
        nul = [r for r in rows if r["cell"] == c and r["dep"] == "null"][0]; q1 = [r for r in rows if r["cell"] == c and r["dep"] == "quadratic" and abs(float(r["rel"]) - 1) < 1e-9][0]
        s = num(nul, t); a = num(q1, t + "_adj")
        if t.startswith("mle") and (mle_ok.get(c) or 0) < .8: vals += ["---", "---"]
        else: vals += [f(s), f(a)]
    il.append("%s & %s \\\\" % (lab, " & ".join(vals)))
inv_tab = ("\\begin{tabular}{lcccccccc}\n\\toprule\n& " + head + " \\\\\n" + "".join("\\cmidrule(lr){%d-%d}" % (2 + 2 * i, 3 + 2 * i) for i in range(4)) +
           "\nTest & size & adj.\\ power & size & adj.\\ power & size & adj.\\ power & size & adj.\\ power \\\\\n\\midrule\n" + "\n".join(il) + "\n\\bottomrule\n\\end{tabular}")
s = io.open(P, encoding="utf-8").read()
if s.count("%% PENDING E4d rows") == 2:
    s = s.replace("%% PENDING E4d rows", power_tab, 1).replace("%% PENDING E4d rows", inv_tab, 1)
else:   # idempotent re-fill: replace the tabular block that follows each label
    for lab, tab in [("tab:power", power_tab), ("tab:invalid", inv_tab)]:
        pat = re.compile(r"(\\label\{" + lab + r"\}\n)\\begin\{tabular\}.*?\\end\{tabular\}", re.S)
        assert len(pat.findall(s)) == 1, lab
        s = pat.sub(lambda m: m.group(1) + tab, s)
s = re.sub(r"(?m)^%% PENDING E4d.*\n", "", s)
io.open(P, "w", encoding="utf-8", newline="\n").write(s)
print("tables filled; PENDING E4d markers removed"); print("\n".join(il))
