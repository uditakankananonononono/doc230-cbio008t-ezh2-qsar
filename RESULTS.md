# CBIO008T results (PREREG.md committed before scoring; src/fetch.py, src/run.py; results/results.json). CPU only, scikit-learn/RDKit. Single seed, no variance estimates.
D1 data: 2559 compounds, 5265 IC50 measurements, 36 patents; mean pIC50 7.3 (sd 1.2); largest patent holds 15% of compounds. 60% of compounds have pIC50 >=7 (potent-skewed, patent data).
D2 random 5-fold, Spearman / RMSE: ridge 0.912/0.489, RF 0.917/0.469, MLP 0.896/0.531. Optimistic, analogue series leak.
D3 leave-patent-out: ridge 0.569/1.041, RF 0.522/1.040, MLP 0.651/1.022. Gate (best >=0.50): met by MLP. Random split overstates by about 0.3 Spearman.
D4 'deep learning helps' gate (MLP beats RF by >=0.03 LPO): met (0.651 vs 0.522), but single seed and the RF is a fixed-hyperparameter baseline; treat as suggestive.
D5 applicability domain (RF, LPO) by max Tanimoto to training: >=0.7 n=562 Spearman 0.777; 0.5-0.7 n=1379 0.692; <0.5 n=618 -0.146 (no predictive ability). Model works only near known chemistry.
D6 potent (pIC50>=7) classification, LPO: AUROC ridge 0.783, RF 0.766, MLP 0.818 (base rate 0.601). Top-10% predicted contain 15-16% of all actives (max possible 16.7%), i.e. about 98% precision at 10% - driven by the high base rate.
D7 fingerprints (RF, LPO Spearman): RDKit topological 0.601, Morgan2 bits 0.522, Morgan3 counts 0.482, MACCS 0.433. Fingerprint choice matters as much as model choice.
D8 noise floor: 993 compounds measured in >=2 assays agree at Spearman 0.969, RMSE 0.248 pIC50. Model LPO RMSE (~1.0) is 4x the assay noise, so the error is model error not label noise.
D9 known drugs: GSK126 is in the data (measured pIC50 8.79); under LPO the RF predicted 6.49 (2.3 log units low). Tazemetostat, EPZ005687 and CPI-1205 are not in this dataset, so unchecked.
D10 scaffold split (RF): Spearman 0.911, RMSE 0.498 - almost the same as random; Bemis-Murcko scaffolds do not stop analogue leakage here, so scaffold split is not an adequate OOD test on this data.
Conclusion: a fingerprint QSAR for EZH2/PRC2 IC50 reaches Spearman about 0.5-0.65 across patent series and has no skill for new chemistry (<0.5 similarity). No prospective compound was designed or validated; nothing here is a drug-discovery result.
Limits: PubChem patent-derived records (assay formats mixed, PRC2 vs EZH2 not curated), median over measurements, one seed, no hyperparameter search, no 3D or docking, no graph neural net. ChEMBL activity endpoint was returning 500 errors at run time, so ChEMBL was not used.
