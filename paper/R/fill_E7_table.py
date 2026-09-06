# Fill Table 4 (tab:penalty) of paper.tex from results/E7_penalty_summary.csv -- item F, any smooth penalty.
# Prints the numbers so the closing sentence of Section 5.6 can be written from them; idempotent.
import csv, io, os, re
HERE = os.path.dirname(os.path.abspath(__file__)); PAPER = os.path.join(HERE, "..", "paper.tex")
CSV = os.path.join(HERE, "..", "..", "results", "E7_penalty_summary.csv")
assert os.path.exists(CSV), "E7_penalty_summary.csv not there yet"
rows = list(csv.DictReader(io.open(CSV, encoding="utf-8"))); R = rows[0]["R"]
pen_lab = {"ridge": "ridge", "gridge": "generalised ridge", "firth": "Firth"}
cell_lab = {"A_p5": "A", "B_p100": r"B, $\kappa = 0.25$"}
def f3(x): return f"{float(x):.3f}"
body = []
for c in ["A_p5", "B_p100"]:
    for p in ["ridge", "gridge", "firth"]:
        r = next(x for x in rows if x["cell"] == c and x["pen"] == p)
        fails = int(r["fails"]); note = "" if fails == 0 else f"$^{{\\dagger}}$"
        body.append(f"{cell_lab[c]} & {pen_lab[p]}{note} & {f3(r['size_unc'])} & {f3(r['size_R0'])} & {f3(r['size_R2'])} & {f3(r['esize_unc'])} & {f3(r['esize_R0'])} & {f3(r['esize_R2'])} & {float(r['defl_R0']):.2f} \\\\")
        print(c, p, "fails", fails, "dec unc/R0/R2", f3(r["size_unc"]), f3(r["size_R0"]), f3(r["size_R2"]), "edge", f3(r["esize_unc"]), f3(r["esize_R0"]), f3(r["esize_R2"]), "defl", f"{float(r['defl_R0']):.3f}")
    if c == "A_p5": body.append(r"\addlinespace")
table = "\n".join([
r"\begin{table}[!t]", r"\centering", r"\small",
r"\caption{Rejection rate of a correctly specified model at nominal level $0.05$ under three smooth penalties, designs A ($\lambda = 100$) " +
r"and B ($\kappa = 0.25$, $\lambda = 416$; the generalised ridge uses $\lambda_j = \lambda(1 + j/p)$): the uncorrected statistic under the " +
r"maximum likelihood covariance, the corrected statistic under the first-order reference $\Om_{\mathrm{MLE}}$, and under the variance at the " +
r"debiased fit $\Sigma_{\mathrm{d}}$; last column, the deflation ratio $\mathrm{E}(\mathrm{SC.HL})/\tr(\Om_{\mathrm{MLE}})$. $" + R + r"$ replications per cell.}",
r"\label{tab:penalty}",
r"\begin{tabular}{llccccccc}", r"\toprule",
r"& & \multicolumn{3}{c}{Decile basis} & \multicolumn{3}{c}{EDGE basis} & \\",
r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}",
r"Design & penalty & uncorrected & first order & $\Sigma_{\mathrm{d}}$ & uncorrected & first order & $\Sigma_{\mathrm{d}}$ & deflation \\", r"\midrule",
*body, r"\bottomrule", r"\end{tabular}", r"\end{table}"])
s = io.open(PAPER, encoding="utf-8").read()
pat = re.compile(r"\\begin\{table\}\[!t\]\n\\centering\n(?:\\small\n)?\\caption\{(?:(?!\\end\{table\}).)*?\\label\{tab:penalty\}.*?\\end\{table\}", re.S)
assert len(pat.findall(s)) == 1, "tab:penalty environment"
s = pat.sub(lambda m: table, s)
s = re.sub(r"%% PENDING E7:.*?\n(?:%%.*?\n)*", "", s)
io.open(PAPER, "w", encoding="utf-8", newline="\n").write(s)
print("Table 4 filled from E7; PENDING E7 removed.")
