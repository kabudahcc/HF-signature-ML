# step39: Surrogate NicheNet ligand prioritization using OmniPath priors
# Receiver: Fibroblast (HF-up genes) and Endothelial (HF-down genes, FCN3 focus)
# Prior: ligand->receptor (OmniPath ligrecextra) + receptor->TF (<=2 steps, OmniPath signaling)
#        + TF->target (CollecTRI). Ligand activity = Fisher enrichment of predicted targets
#        among receiver HF-DE genes vs expressed background.
import json, csv, math, os
from collections import defaultdict, deque

BASE = "D:/早早生信分析学习课程/新R文件"
PROJ = BASE + "/for_publication/Paper2_HF_signature_ML"
MODEL = PROJ + "/data_external/nch/model"
DEDIR = BASE + "/results_scRNA/step2_DEG_composition"
OUT = PROJ + "/revision_round2"

try:
    from scipy.stats import fisher_exact
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False
    def logcomb(n, k):
        if k < 0 or k > n: return float("-inf")
        return math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1)
    def fisher_exact(table, alternative="greater"):
        # one-sided greater (enrichment) hypergeometric p-value
        (a, b), (c, d) = table
        n = a + b + c + d; r1 = a + b; c1 = a + c
        lo = max(0, r1 - (n - c1)); hi = min(r1, c1)
        denom = logcomb(n, r1)
        p = sum(math.exp(logcomb(c1, x) + logcomb(n - c1, r1 - x) - denom)
                for x in range(a, hi + 1) if lo <= x <= hi)
        return None, min(p, 1.0)

def load_json(name):
    with open(os.path.join(MODEL, name)) as f:
        return json.load(f)

# ---------- 1. Ligand-receptor pairs ----------
lr_raw = load_json("omnipath_ligrecextra_raw.json")
lr_pairs = defaultdict(lambda: {"n_sources": 0})
for e in lr_raw:
    L, R = e["source_genesymbol"], e["target_genesymbol"]
    if L and R and L != R:
        key = (L, R)
        lr_pairs[key]["n_sources"] = max(lr_pairs[key]["n_sources"], len(e.get("sources", [])))
receptors_of = defaultdict(set)   # ligand -> set(receptors)
for (L, R) in lr_pairs:
    receptors_of[L].add(R)
print(f"LR pairs: {len(lr_pairs)}, ligands: {len(receptors_of)}")

# ---------- 2. Signaling graph (receptor -> TF, <=2 steps) ----------
sig_edges = defaultdict(set)
with open(os.path.join(MODEL, "omnipath_signaling.tsv")) as f:
    r = csv.DictReader(f, delimiter="\t")
    for e in r:
        if e["is_directed"] == "True":
            s, t = e["source_genesymbol"], e["target_genesymbol"]
            if s and t and not s.startswith("COMPLEX:") and not t.startswith("COMPLEX:"):
                sig_edges[s].add(t)
print(f"signaling edges: {sum(len(v) for v in sig_edges.values())}, nodes: {len(sig_edges)}")

# ---------- 3. CollecTRI TF -> target ----------
tf_targets = defaultdict(set)
with open(os.path.join(MODEL, "collectri.tsv")) as f:
    r = csv.DictReader(f, delimiter="\t")
    for e in r:
        tf, tg = e["source_genesymbol"], e["target_genesymbol"]
        if tf and tg:
            tf_targets[tf].add(tg)
TF_SET = set(tf_targets)
print(f"TFs: {len(TF_SET)}, TF-target edges: {sum(len(v) for v in tf_targets.values())}")

def reachable_tfs(receptor, max_depth=2):
    if receptor in TF_SET:
        return {receptor}
    seen, frontier = {receptor}, {receptor}
    for _ in range(max_depth):
        nxt = set()
        for node in frontier:
            nxt |= sig_edges.get(node, set())
        nxt -= seen
        tfs = nxt & TF_SET
        if tfs:
            return tfs
        seen |= nxt
        frontier = nxt
    return set()

# ---------- 4. DE tables ----------
CELL_TYPES = ["Fibroblast", "Endothelial", "Macrophage", "CM", "Pericyte_SMC",
              "T_cell", "NK", "B_cell", "Mast", "Neutrophil"]
# filename mapping
def de_path(ct, comp):
    return os.path.join(DEDIR, f"DE_{ct}_{comp}.csv")

de = {}  # (ct, comp) -> {gene: dict(log2fc, p_adj, pct1, pct2)}
for ct in CELL_TYPES:
    for comp in ["DCM_vs_N", "ICM_vs_N", "DCM_vs_ICM"]:
        p = de_path(ct, comp)
        if not os.path.exists(p):
            continue
        tab = {}
        with open(p) as f:
            for row in csv.DictReader(f):
                try:
                    tab[row["gene"]] = {
                        "log2fc": float(row["avg_log2FC"]),
                        "p_adj": float(row["p_val_adj"]),
                        "pct1": float(row["pct.1"]),
                        "pct2": float(row["pct.2"]),
                    }
                except (ValueError, KeyError):
                    pass
        de[(ct, comp)] = tab
print(f"DE tables loaded: {len(de)}")

def expr_frac(ct, gene):
    """max expression fraction across conditions for a gene in a cell type"""
    best = 0.0
    for comp in ["DCM_vs_N", "ICM_vs_N", "DCM_vs_ICM"]:
        tab = de.get((ct, comp))
        if tab and gene in tab:
            best = max(best, tab[gene]["pct1"], tab[gene]["pct2"])
    return best

def hf_de_set(ct, direction):
    """HF (DCM or ICM) vs N DE genes: direction 'up' or 'down'"""
    genes = set()
    for comp in ["DCM_vs_N", "ICM_vs_N"]:
        tab = de.get((ct, comp), {})
        for g, s in tab.items():
            if s["p_adj"] < 0.05 and abs(s["log2fc"]) > 0.25:
                if (direction == "up" and s["log2fc"] > 0) or (direction == "down" and s["log2fc"] < 0):
                    genes.add(g)
    return genes

def background_genes(ct):
    bg = set()
    for comp in ["DCM_vs_N", "ICM_vs_N", "DCM_vs_ICM"]:
        bg |= set(de.get((ct, comp), {}).keys())
    return bg

# ---------- 5. Ligand activity (Fisher) ----------
def ligand_activity(receiver_ct, direction, ligand_cond="any"):
    geneset = hf_de_set(receiver_ct, direction) & background_genes(receiver_ct)
    bg = background_genes(receiver_ct)
    print(f"[{receiver_ct}/{direction}] geneset={len(geneset)}, background={len(bg)}")
    rows = []
    for L, recs in receptors_of.items():
        # receptor expressed in receiver (pct>=0.1 any condition)
        recs_expr = [R for R in recs if expr_frac(receiver_ct, R) >= 0.1]
        if not recs_expr:
            continue
        # ligand expressed in >=1 sender cell type (any condition)
        senders = [ct for ct in CELL_TYPES if ct != receiver_ct and expr_frac(ct, L) >= 0.1]
        auto = expr_frac(receiver_ct, L) >= 0.1
        if not senders and not auto:
            continue
        # predicted targets
        tfs = set()
        for R in recs_expr:
            tfs |= reachable_tfs(R)
        targets = set()
        for tf in tfs:
            targets |= tf_targets[tf]
        targets &= bg
        if len(targets) < 5:
            continue
        a = len(targets & geneset)
        b = len(targets) - a
        c = len(geneset) - a
        d = len(bg) - len(geneset) - b
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({"ligand": L, "receptors": ";".join(sorted(recs_expr)),
                     "n_targets_bg": len(targets), "n_targets_in_geneset": a,
                     "frac_in_geneset": round(a / max(len(geneset), 1), 4),
                     "p_fisher": p, "senders": ";".join(senders),
                     "autocrine": auto, "n_tfs": len(tfs)})
    # BH adjust
    rows.sort(key=lambda r: r["p_fisher"])
    m = len(rows)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, rows[i]["p_fisher"] * m / (i + 1))
        rows[i]["p_adj_BH"] = min(prev, 1.0)
    return rows, geneset, bg

def save_rows(rows, path, extra_cols=()):
    cols = ["ligand", "receptors", "n_targets_bg", "n_targets_in_geneset",
            "frac_in_geneset", "p_fisher", "p_adj_BH", "senders", "autocrine", "n_tfs"] + list(extra_cols)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})

# ---------- 6. Run: Fibroblast receiver (HF-up) ----------
fib_rows, fib_geneset, fib_bg = ligand_activity("Fibroblast", "up")
save_rows(fib_rows, os.path.join(OUT, "niche_surrogate_fibroblast_ligand_activity.csv"))
print("top fibroblast ligands:")
for r in fib_rows[:15]:
    print(f"  {r['ligand']:10s} rec={r['receptors'][:40]:40s} hits={r['n_targets_in_geneset']:3d}/{r['n_targets_bg']:4d} p={r['p_fisher']:.2e} padj={r['p_adj_BH']:.2e}")

# check our signature genes among targets of top ligands
SIG = ["ASPN", "ISLR", "HMGN2", "PLTP"]
print("signature genes in fibroblast HF-up set:", [g for g in SIG if g in fib_geneset])

# ---------- 7. Sender evidence for top fibroblast ligands ----------
evidence = []
for r in fib_rows[:20]:
    L = r["ligand"]
    for ct in CELL_TYPES:
        if ct == "Fibroblast":
            continue
        for comp in ["DCM_vs_N", "ICM_vs_N"]:
            tab = de.get((ct, comp), {})
            if L in tab:
                s = tab[L]
                evidence.append({"ligand": L, "sender": ct, "comparison": comp,
                                 "log2fc": s["log2fc"], "p_adj": s["p_adj"],
                                 "pct_hf": s["pct1"], "pct_n": s["pct2"]})
with open(os.path.join(OUT, "niche_surrogate_fibroblast_top20_sender_evidence.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["ligand", "sender", "comparison", "log2fc", "p_adj", "pct_hf", "pct_n"])
    w.writeheader(); w.writerows(evidence)

# ---------- 8. Endothelial receiver (HF-down; FCN3) ----------
endo_rows, endo_geneset, _ = ligand_activity("Endothelial", "down")
save_rows(endo_rows, os.path.join(OUT, "niche_surrogate_endothelial_ligand_activity.csv"))
print("FCN3 in endothelial HF-down set:", "FCN3" in endo_geneset)
print("top endothelial (down) ligands:")
for r in endo_rows[:10]:
    print(f"  {r['ligand']:10s} hits={r['n_targets_in_geneset']:3d}/{r['n_targets_bg']:4d} p={r['p_fisher']:.2e} padj={r['p_adj_BH']:.2e}")

# which prioritized ligands have FCN3 as predicted target?
fcn3_ligands = []
for L, recs in receptors_of.items():
    recs_expr = [R for R in recs if expr_frac("Endothelial", R) >= 0.1]
    if not recs_expr:
        continue
    tfs = set()
    for R in recs_expr:
        tfs |= reachable_tfs(R)
    targets = set()
    for tf in tfs:
        targets |= tf_targets[tf]
    if "FCN3" in targets:
        # is ligand down in HF in any sender?
        down_ev = []
        for ct in CELL_TYPES:
            for comp in ["DCM_vs_N", "ICM_vs_N"]:
                tab = de.get((ct, comp), {})
                if L in tab and tab[L]["p_adj"] < 0.05 and tab[L]["log2fc"] < -0.25:
                    down_ev.append(f"{ct}/{comp}:{tab[L]['log2fc']:.2f}")
        fcn3_ligands.append({"ligand": L, "receptors": ";".join(sorted(recs_expr)),
                             "down_in_HF": ";".join(down_ev)})
with open(os.path.join(OUT, "niche_surrogate_FCN3_upstream_ligands.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["ligand", "receptors", "down_in_HF"])
    w.writeheader(); w.writerows(fcn3_ligands)
print(f"ligands with FCN3 as predicted target: {len(fcn3_ligands)}, of which down in HF: {sum(1 for x in fcn3_ligands if x['down_in_HF'])}")

# ---------- 9. summary ----------
with open(os.path.join(OUT, "niche_surrogate_summary.txt"), "w") as f:
    f.write("Surrogate NicheNet (OmniPath priors) — summary\n")
    f.write(f"LR pairs: {len(lr_pairs)}; ligands: {len(receptors_of)}; TFs: {len(TF_SET)}\n")
    f.write(f"Fibroblast HF-up geneset: {len(fib_geneset)} (bg {len(fib_bg)})\n")
    f.write(f"Endothelial HF-down geneset: {len(endo_geneset)}\n")
    f.write("Signature genes in fibroblast HF-up: " + ",".join(g for g in SIG if g in fib_geneset) + "\n")
    f.write("Top-10 fibroblast ligands:\n")
    for r in fib_rows[:10]:
        f.write(f"  {r['ligand']}\t{r['receptors']}\t{r['n_targets_in_geneset']}/{r['n_targets_bg']}\tp={r['p_fisher']:.3e}\tBH={r['p_adj_BH']:.3e}\n")
    f.write(f"FCN3 predicted upstream ligands: {len(fcn3_ligands)} (down in HF: {sum(1 for x in fcn3_ligands if x['down_in_HF'])})\n")
print("DONE. scipy:", HAVE_SCIPY)
