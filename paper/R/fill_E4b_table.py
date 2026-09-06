# Fill Table 3 (tab:power) of paper.tex from results/E4b_final_summary.csv -- the definitive CALM-vs-prepivot grid.
# Run once E4b has written its summary; idempotent (replaces the whole table environment + removes the PENDING comment).
import csv, io, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PAPER = os.path.join(HERE, "..", "paper.tex")
CSV = os.path.join(HERE, "..", "..", "results", "E4b_final_summary.csv")
assert os.path.exists(CSV), "E4b_final_summary.csv not there yet"
rows = list(csv.DictReader(io.open(CSV, encoding="utf-8")))
R, B = rows[0]["R"], rows[0]["B"]
cell_lab = {"A_p5": r"A ($\kappa = 0.01$)", "B_p40": r"B, $\kappa = 0.10$", "B_p100": r"B, $\kappa = 0.25$", "B_p160": r"B, $\kappa = 0.40$"}
dep_lab = {"null": "---", "index": "cubic in the index", "offindex": r"$x_1x_2$ off the index"}
order = ["A_p5", "B_p40", "B_p100", "B_p160"]; dorder = [("null", "0"), ("index", "0.4"), ("index", "0.7"), ("offindex", "0.4"), ("offindex", "0.7")]
def f3(x): return f"{float(x):.3f}"
def diff(s):  # "+0.012 (0.005)" -> "$+0.012$ (0.005)"
    m = re.match(r"([+-]?[0-9.]+) \(([0-9.]+)\)", s); return f"${m.group(1)}$ ({m.group(2)})"
body = []
for c in order:
    first = True
    for dep, rel in dorder:
        r = next(x for x in rows if x["cell"] == c and x["dep"] == dep and float(x["rel"]) == float(rel))
        lab = cell_lab[c] if first else ""; first = False
        body.append(f"{lab} & {dep_lab[dep]} & {rel if dep != 'null' else '0'} & {r['deg']} & {f3(r['cf_hl'])} & {f3(r['pp_hl'])} & {diff(r['diff_hl'])} & {f3(r['cf_ad'])} & {f3(r['pp_edge'])} & {diff(r['diff_ad'])} \\\\")
    body.append(r"\addlinespace")
body = body[:-1]
table = "\n".join([
r"\begin{table}[!t]", r"\centering", r"\small",
r"\caption{Size and power at nominal level $0.05$ of CALM and of prepivoting with $" + B + r"$ refits, paired on identical datasets, in the four " +
r"$\kappa$ cells of Section~\ref{sec:designs}. Departures are a Hermite cubic along the true index and an omitted interaction $x_1x_2$ off the " +
r"index, each at relative size $0.4$ and $0.7$ of the linear predictor; the rows with relative size $0$ are size. $k$ is the degree chosen by the " +
r"adaptive rule in the majority of replications. Diff is the paired difference CALM minus prepivot with its standard error. $" + R + r"$ replications per cell.}",
r"\label{tab:power}",
r"\begin{tabular}{lllccccccc}", r"\toprule",
r"& & & & \multicolumn{3}{c}{Decile basis (SC.HL)} & \multicolumn{3}{c}{Adaptive EDGE basis} \\",
r"\cmidrule(lr){5-7}\cmidrule(lr){8-10}",
r"Design & departure & size & $k$ & CALM & prepivot & diff (s.e.) & CALM & prepivot & diff (s.e.) \\", r"\midrule",
*body, r"\bottomrule", r"\end{tabular}", r"\end{table}"])
s = io.open(PAPER, encoding="utf-8").read()
pat = re.compile(r"\\begin\{table\}\[!t\]\n\\centering\n(?:\\small\n)?\\caption\{(?:(?!\\end\{table\}).)*?\\label\{tab:power\}.*?\\end\{table\}", re.S)
assert len(pat.findall(s)) == 1, "tab:power environment"
s = pat.sub(lambda m: table, s)
s = re.sub(r"%% PENDING E4b:.*?\n(?:%%.*?\n)*", "", s)
io.open(PAPER, "w", encoding="utf-8", newline="\n").write(s)
print("Table 3 filled from E4b; PENDING E4b removed. Rows:"); print("\n".join(body))
