# Series B (CALM) manuscript audit -- RSS house style + build health.  Prints the FULL report (never head it).
import io, os, re, subprocess, sys, glob
P = r"C:\Users\ebrah\.gemini\Projects\PDFs\paper_seriesB_closedform\paper"
os.chdir(P)
tex = io.open("paper.tex", encoding="utf-8").read()
body = re.sub(r"(?m)^%.*$", "", tex)              # drop comment lines
fails = []
def check(ok, msg):
    print(("PASS  " if ok else "FAIL  ") + msg)
    if not ok: fails.append(msg)

# 1. build
for _ in range(2):
    subprocess.run(["pdflatex", "-interaction=nonstopmode", "paper.tex"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
log = io.open("paper.log", encoding="latin-1").read()
errs = len(re.findall(r"(?m)^!", log)); undef = len(re.findall(r"(Citation|Reference) .* undefined", log))
pages = int(re.findall(r"\((\d+) pages", log)[-1]); over = log.count("Overfull")
check(errs == 0, f"pdflatex errors = {errs}"); check(undef == 0, f"undefined cites/refs = {undef}")
check(pages <= 24, f"pages = {pages} (RSS: >24 unlikely to be accepted)"); check(over == 0, f"overfull boxes = {over}")

# 2. abstract
ab = re.search(r"\\abstract\{(.*?)\}\s*\n\s*\n", body, re.S).group(1)
ab_words = len(re.sub(r"\$[^$]*\$", "X", ab).split())
check(ab_words <= 200, f"abstract words = {ab_words}")
check("\\cite" not in ab, "abstract has no citations")
check(not re.search(r"\b[A-Z]{2,}\b(?<!CALM)", ab.replace("CALM", "")), "abstract has no abbreviations (CALM is the name)")

# 3. keywords alphabetical + each in abstract
kw = [k.strip() for k in re.search(r"\\keywords\{(.*?)\}", body, re.S).group(1).replace("\n", " ").split(";")]
check(kw == sorted(kw, key=str.lower), f"keywords alphabetical: {kw}")
check(5 <= len(kw) <= 6, f"keyword count = {len(kw)} (5-6)")
abl = re.sub(r"\s+", " ", ab.lower()).replace("-", " ")
for k in kw:
    stem = re.sub(r"\s+", " ", k.lower()).replace("-", " ")
    check(stem in abl, f"keyword in abstract: '{k}'")

# 4. house style
check("\\begin{itemize}" not in body and "\\begin{enumerate}" not in body, "no bullet/enumerate lists")
check("\\footnote" not in body, "no footnotes")
check(not re.search(r"\\newtheorem\{\w+\}\[theorem\]", body), "theorems numbered individually by type (no [theorem] counter sharing)")
check("\\newcommand{\\societylogo}{}" in tex, "societylogo defined")
check("\\authormark" not in body, "no \\authormark")
nfig = len(re.findall(r"\\begin\{figure", body)); nalt = body.count("\\figalttext")
check(nfig == nalt and nfig > 0, f"figures = {nfig}, figalttext = {nalt}")
ncap = len(re.findall(r"\\caption\{", body))
print(f"INFO  captions = {ncap} (figures + tables)")
check(not re.search(r"\\includegraphics[^{]*\{[^}]*\.(jpe?g|png)\}", body), "figures are vector (no jpg/png)")
check("\\section*{Data availability}" in body or "Data availability" in body, "Data availability section present")
check(re.search(r"AI|artificial intelligence|language model", body, re.I) is not None, "AI-use disclosure present")
check("\\hline" not in body, "no \\hline (booktabs only, no vertical/inner rules)")
check("|" not in "".join(re.findall(r"\\begin\{tabular\}\{([^}]*)\}", body)), "no vertical table lines")
# UK spelling: RSS requires it
main_body = body.split("\\begin{thebibliography}")[0]
main_body = re.sub(r"\\(includegraphics|label|ref)\{[^}]*\}", "", main_body)   # file names and labels are not prose
us = sorted(set(re.findall(r"\b\w*(?:ize|izes|ized|izing|ization|izations|behavior|center|centered|modeling|modeled|color|colors|favor|analyze|analyzed|analyzes|analyzing)\b", main_body)))
us = [w for w in us if not w.lower().endswith(("size", "sizes", "seize", "prize", "prizes"))]
check(not us, f"US spellings (RSS wants UK): {us}")
# pending markers
pend = re.findall(r"(?m)^%% PENDING.*$", tex)
print(f"INFO  PENDING markers ({len(pend)}): {pend}")
# running head
rh = re.search(r"\\title\[(.*?)\]", body, re.S).group(1)
check(len(rh) <= 50, f"running head length = {len(rh)} (<=50): '{rh}'")
# bibitems all cited, all cites have bibitems
keys = set(re.findall(r"\\bibitem\[[^\]]*\]\{([^}]+)\}", body))
cited = set(); [cited.update(c.strip() for c in m.split(",")) for m in re.findall(r"\\cite(?:t|p|alp|alt|author|year)?\*?(?:\[[^\]]*\])*\{([^}]+)\}", body)]
check(cited <= keys, f"cites without bibitem: {sorted(cited - keys)}")
check(keys <= cited, f"bibitems never cited: {sorted(keys - cited)}")
print(f"INFO  bibitems = {len(keys)}, distinct cites = {len(cited)}")
# tables/figures referenced
for lab in re.findall(r"\\label\{((?:tab|fig):[^}]+)\}", body):
    check(f"\\ref{{{lab}}}" in body, f"{lab} referenced in text")
# theorem sanity in PDF text: Proposition 1 must exist if propositions used
if "\\begin{proposition}" in body:
    txt = subprocess.run(["pdftotext", "paper.pdf", "-"], capture_output=True, text=True, encoding="utf-8", errors="ignore").stdout
    check("Proposition 1" in txt, "PDF shows 'Proposition 1' (own counter)")
# results files the paper depends on
for f in ["E2b_deflation_summary.csv", "E4b_final_summary.csv", "E7_penalty_summary.csv", "E9_kappa_inflation_summary.csv"]:
    hit = glob.glob(os.path.join("..", "results", f.replace("_summary", "*")))
    print(f"INFO  results/{f}: {'present' if hit else 'MISSING'}")
print(f"\n=== {len(fails)} FAIL(s) ===")
for f in fails: print("  -", f)
sys.exit(1 if fails else 0)
