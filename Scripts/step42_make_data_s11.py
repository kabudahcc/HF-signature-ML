# build Supplementary_Data_S11.xlsx — surrogate ligand prioritization tables
import pandas as pd, os

PROJ = "D:/早早生信分析学习课程/新R文件/for_publication/Paper2_HF_signature_ML"
OUT = PROJ + "/revision_round2"
DEST = PROJ + "/Supplementary/Supplementary_Data_S11.xlsx"

fib = pd.read_csv(OUT + "/niche_surrogate_fibroblast_ligand_activity_hiconf.csv")
ev = pd.read_csv(OUT + "/niche_surrogate_fibroblast_top20_sender_evidence.csv")
endo = pd.read_csv(OUT + "/niche_surrogate_endothelial_ligand_activity.csv")
notes = pd.DataFrame({"note": [
    "Surrogate NicheNet-style ligand prioritization using OmniPath prior knowledge (the NicheNet ligand-target prior model, Zenodo record 7074291, was unavailable in our computing environment).",
    "Prior graph: ligand->receptor (OmniPath ligand-receptor extra resource; pairs supported by >=2 resources), receptor->TF (directed OmniPath signaling network, <=2 intermediates), TF->target (CollecTRI regulons).",
    "Receiver: fibroblast HF-upregulated genes (DCM vs. N or ICM vs. N, |log2FC|>0.25, FDR<0.05; n=441) vs. fibroblast-expressed background (n=2,184); one-sided Fisher exact test, Benjamini-Hochberg FDR.",
    "Sheet 1: full ranked fibroblast-receiver ligand activity (high-confidence ligands only). Sheet 2: sender-cell differential expression of the top-20 ligands. Sheet 3: endothelial-receiver ligand activity against HF-DOWN genes (exploratory).",
    "Intracellular proteins mis-annotated as ligands in ligand-receptor resources (VIM, PKM, TYROBP) appear in the ranked list but were excluded from interpretation.",
    "FCN3 is absent from CollecTRI regulons; no expressed ligand could be linked to FCN3 regulation, so upstream control of FCN3 remains unresolved.",
    "Sources: OmniPath (omnipathdb.org, interactions downloaded 2026-08-05; Turei et al. Mol Syst Biol 2021); CollecTRI (via OmniPath; Muller-Dott et al. Nucleic Acids Res 2023); NicheNet framework (Browaeys et al. Nat Methods 2020).",
]})

with pd.ExcelWriter(DEST, engine="openpyxl") as w:
    fib.to_excel(w, sheet_name="fibroblast_ligand_activity", index=False)
    ev.to_excel(w, sheet_name="sender_evidence_top20", index=False)
    endo.to_excel(w, sheet_name="endothelial_down_activity", index=False)
    notes.to_excel(w, sheet_name="notes", index=False)
print("written", DEST, os.path.getsize(DEST), "bytes;", len(fib), len(ev), len(endo), "rows")
