# step40: refine surrogate NicheNet results
# - confidence filter (LR pairs supported by >=2 resources)
# - APOE axis detail + signature-gene target coverage
# - FCN3 in CollecTRI check
# - Figure S23 (fibroblast ligand prioritization)
import csv, os, json
from collections import defaultdict

BASE = "D:/早早生信分析学习课程/新R文件"
PROJ = BASE + "/for_publication/Paper2_HF_signature_ML"
MODEL = PROJ + "/data_external/nch/model"
OUT = PROJ + "/revision_round2"

# LR pairs with source counts
lr_raw = json.load(open(os.path.join(MODEL, "omnipath_ligrecextra_raw.json")))
pair_src = defaultdict(set)
for e in lr_raw:
    L, R = e["source_genesymbol"], e["target_genesymbol"]
    if L and R and L != R:
        pair_src[(L, R)] |= set(e.get("sources", []))
conf_pairs = {p for p, s in pair_src.items() if len(s) >= 2}
print(f"LR pairs total {len(pair_src)}, high-confidence (>=2 resources) {len(conf_pairs)}")
conf_ligands = sorted({L for L, R in conf_pairs})
print("high-confidence ligands:", len(conf_ligands))

# re-rank activity restricted to high-confidence ligands
act = list(csv.DictReader(open(os.path.join(OUT, "niche_surrogate_fibroblast_ligand_activity.csv"))))
act_conf = [r for r in act if r["ligand"] in set(conf_ligands)]
print(f"activity rows: {len(act)} -> high-confidence {len(act_conf)}")
print("top-12 high-confidence fibroblast ligands:")
for r in act_conf[:12]:
    print(f"  {r['ligand']:8s} rec={r['receptors'][:35]:35s} hits={r['n_targets_in_geneset']}/{r['n_targets_bg']} BH={float(r['p_adj_BH']):.2e}")

# APOE stats
for r in act:
    if r["ligand"] == "APOE":
        print("APOE row:", {k: r[k] for k in ["receptors", "n_targets_bg", "n_targets_in_geneset", "p_fisher", "p_adj_BH", "senders"]})

# ---- signature gene coverage by top axes ----
# rebuild targets for key ligands to check ASPN/ISLR/HMGN2/PLTP
sig_edges = defaultdict(set)
with open(os.path.join(MODEL, "omnipath_signaling.tsv")) as f:
    for e in csv.DictReader(f, delimiter="\t"):
        if e["is_directed"] == "True":
            s, t = e["source_genesymbol"], e["target_genesymbol"]
            if s and t and not s.startswith("COMPLEX:") and not t.startswith("COMPLEX:"):
                sig_edges[s].add(t)
tf_targets = defaultdict(set)
with open(os.path.join(MODEL, "collectri.tsv")) as f:
    for e in csv.DictReader(f, delimiter="\t"):
        tf, tg = e["source_genesymbol"], e["target_genesymbol"]
        if tf and tg:
            tf_targets[tf].add(tg)
TF_SET = set(tf_targets)

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
        seen |= nxt; frontier = nxt
    return set()

SIG = ["ASPN", "ISLR", "HMGN2", "PLTP"]
lig_recs = defaultdict(set)
for (L, R) in conf_pairs:
    lig_recs[L].add(R)
cov_rows = []
for L in [r["ligand"] for r in act_conf[:15]]:
    tfs = set()
    for R in lig_recs[L]:
        tfs |= reachable_tfs(R)
    tgts = set()
    for tf in tfs:
        tgts |= tf_targets[tf]
    hit = [g for g in SIG if g in tgts]
    cov_rows.append({"ligand": L, "signature_targets": ";".join(hit)})
    print(f"  {L:8s} predicts signature genes: {hit}")
with open(os.path.join(OUT, "niche_surrogate_signature_target_coverage.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["ligand", "signature_targets"]); w.writeheader(); w.writerows(cov_rows)

# ---- FCN3 in CollecTRI ----
fcn3_tfs = [tf for tf, tg in tf_targets.items() if "FCN3" in tg]
print("TFs targeting FCN3 in CollecTRI:", fcn3_tfs if fcn3_tfs else "NONE (FCN3 absent from prior)")

# save high-confidence activity table
cols = act[0].keys()
with open(os.path.join(OUT, "niche_surrogate_fibroblast_ligand_activity_hiconf.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(act_conf)
print("saved hiconf table")
