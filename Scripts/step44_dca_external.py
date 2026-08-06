# step44: decision curve analysis for the two largest scale-compatible external cohorts
import csv, os, sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
import numpy as np
import matplotlib.pyplot as plt

PROJ = "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
rows = list(csv.DictReader(open(PROJ + "/Data/HF_signature_external_validation_predictions.csv")))
COHORTS = ["GSE116250", "GSE55296"]

def dca(prob, y, thresholds):
    n = len(y); prev = np.mean(y)
    nb_model, nb_all = [], []
    for pt in thresholds:
        pred = prob >= pt
        tp = np.sum(pred & (y == 1)); fp = np.sum(pred & (y == 0))
        nb_model.append(tp / n - fp / n * (pt / (1 - pt)))
        nb_all.append(prev - (1 - prev) * (pt / (1 - pt)))
    return np.array(nb_model), np.array(nb_all)

thresholds = np.linspace(0.02, 0.80, 200)
setup_plot()
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
out_rows = []
for ax, co in zip(axes, COHORTS):
    sub = [r for r in rows if r["Dataset"] == co]
    prob = np.array([float(r["Prob_Disease"]) for r in sub])
    y = np.array([1 if r["TrueLabel"] == "Disease" else 0 for r in sub])
    nb_m, nb_a = dca(prob, y, thresholds)
    for pt, a, b in zip(thresholds, nb_m, nb_a):
        out_rows.append({"cohort": co, "threshold": round(pt, 3),
                         "net_benefit_model": round(a, 4), "net_benefit_treat_all": round(b, 4)})
    ax.plot(thresholds, nb_m, color="#b2182b", lw=2, label="39-gene model")
    ax.plot(thresholds, nb_a, color="grey", lw=1.2, ls="--", label="Treat all")
    ax.axhline(0, color="black", lw=1.2, ls=":", label="Treat none")
    ax.set_title(f"{co}  (n = {len(y)}; {int(y.sum())} HF / {int((1-y).sum())} non-failing)", fontsize=10)
    ax.set_xlabel("Threshold probability")
    ax.set_ylim(-0.15, 0.65)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("Net benefit")
fig.suptitle("Decision curve analysis of the frozen model in external cohorts", y=1.02, fontsize=11)
fig.savefig(PROJ + "/revision_round2/S24_dca_external.png", dpi=300, bbox_inches="tight")

with open(PROJ + "/revision_round2/dca_external_curves.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["cohort", "threshold", "net_benefit_model", "net_benefit_treat_all"])
    w.writeheader(); w.writerows(out_rows)

# summarize the threshold range where model beats both strategies
for co in COHORTS:
    sub = [r for r in rows if r["Dataset"] == co]
    prob = np.array([float(r["Prob_Disease"]) for r in sub])
    y = np.array([1 if r["TrueLabel"] == "Disease" else 0 for r in sub])
    nb_m, nb_a = dca(prob, y, thresholds)
    win = (nb_m > nb_a) & (nb_m > 0)
    if win.any():
        ts = thresholds[win]
        print(f"{co}: model superior over thresholds {ts.min():.2f}-{ts.max():.2f}")
    else:
        print(f"{co}: no threshold range with superiority")
print("saved S24_dca_external.png + dca_external_curves.csv")
