# step41: Figure S23 — surrogate ligand prioritization for fibroblast HF program
# Panel A: ligand activity (-log10 BH p, Fisher enrichment of predicted targets in fibroblast HF-up genes)
# Panel B: sender-cell differential expression of ligands (HF vs N, max signed log2FC across DCM/ICM)
import csv, os, sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
import matplotlib.pyplot as plt
import numpy as np

BASE = "D:/早早生信分析学习课程/新R文件"
PROJ = BASE + "/for_publication/Paper2_HF_signature_ML"
OUT = PROJ + "/revision_round2"
DEDIR = BASE + "/results_scRNA/step2_DEG_composition"

EXCLUDE = {"VIM", "PKM", "TYROBP"}  # intracellular proteins mis-annotated as ligands in LR resources

act = list(csv.DictReader(open(os.path.join(OUT, "niche_surrogate_fibroblast_ligand_activity_hiconf.csv"))))
act = [r for r in act if r["ligand"] not in EXCLUDE]
# ensure APOE included
top = act[:11]
ligs = [r["ligand"] for r in top]
if "APOE" not in ligs:
    apoe = [r for r in act if r["ligand"] == "APOE"][0]
    top.append(apoe); ligs.append("APOE")
rec_map = {r["ligand"]: r["receptors"] for r in top}
bh_map = {r["ligand"]: float(r["p_adj_BH"]) for r in top}
hit_map = {r["ligand"]: f"{r['n_targets_in_geneset']}/{r['n_targets_bg']}" for r in top}

CELLS = ["Macrophage", "Endothelial", "CM", "Pericyte_SMC", "Fibroblast", "T_cell", "NK", "B_cell", "Neutrophil", "Mast"]
def sender_log2fc(L, ct):
    best = None
    for comp in ["DCM_vs_N", "ICM_vs_N"]:
        p = os.path.join(DEDIR, f"DE_{ct}_{comp}.csv")
        if not os.path.exists(p):
            continue
        for row in csv.DictReader(open(p)):
            if row["gene"] == L:
                try:
                    lfc, pa = float(row["avg_log2FC"]), float(row["p_val_adj"])
                except ValueError:
                    continue
                if pa < 0.05 and (best is None or abs(lfc) > abs(best)):
                    best = lfc
    return best if best is not None else 0.0

M = np.array([[sender_log2fc(L, ct) for ct in CELLS] for L in ligs])

setup_plot()
fig = plt.figure(figsize=(11, 6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 2.2], wspace=0.45)

# Panel A
axA = fig.add_subplot(gs[0])
y = np.arange(len(ligs))[::-1]
vals = [-np.log10(max(bh_map[L], 1e-300)) for L in ligs]
colors = ["#b2182b" if L == "APOE" else "#2166ac" for L in ligs]
axA.barh(y, vals, color=colors, height=0.7)
axA.set_yticks(y)
axA.set_yticklabels([f"{L}  ({rec_map[L].replace(';', '·')})" for L in ligs], fontsize=8)
axA.set_xlabel("Ligand activity  (-log$_{10}$ BH $P$)")
for yi, L in zip(y, ligs):
    axA.text(-np.log10(bh_map[L]) + 0.15, yi, hit_map[L], va="center", fontsize=7, color="0.3")
axA.set_title("A  Prioritized ligands (OmniPath surrogate)", loc="left", fontsize=10)
axA.spines[["top", "right"]].set_visible(False)

# Panel B
axB = fig.add_subplot(gs[1])
vmax = 2.0
im = axB.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
axB.set_xticks(range(len(CELLS)))
axB.set_xticklabels(CELLS, rotation=45, ha="right", fontsize=8)
axB.set_yticks(y[::-1] if False else range(len(ligs)))
axB.set_yticklabels(ligs, fontsize=8)
for i in range(len(ligs)):
    for j in range(len(CELLS)):
        v = M[i, j]
        if abs(v) > 0.01:
            axB.text(j, i, f"{v:+.1f}", ha="center", va="center",
                     fontsize=6.5, color="white" if abs(v) > 1.0 else "0.2")
axB.set_title("B  Ligand expression change in sender cells (HF vs non-failing)", loc="left", fontsize=10)
cb = fig.colorbar(im, ax=axB, fraction=0.03, pad=0.02)
cb.set_label("max signed log$_2$FC (FDR<0.05)", fontsize=8)

fig.savefig(os.path.join(OUT, "S23_niche_surrogate.png"), dpi=300, bbox_inches="tight")
print("saved S23_niche_surrogate.png")

# underlying table for SI
with open(os.path.join(OUT, "niche_surrogate_S23_source_data.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["ligand", "receptors", "hits", "p_adj_BH"] + [f"log2fc_{c}" for c in CELLS])
    for i, L in enumerate(ligs):
        w.writerow([L, rec_map[L], hit_map[L], f"{bh_map[L]:.3g}"] + [f"{v:.2f}" for v in M[i]])
print("saved source data")
