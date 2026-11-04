# GEO myotube signatures

Generated from GEO-deposited human muscle-cell data for the Cell2Sentence analysis. CSVs use the required `gene`, `score`, and `direction` columns. Scores are ranked test statistics, not log2 fold changes; keep `SIGNATURE_SCORE_MODE = "rank_proxy"`.

| GEO study / treatment | Design and exported genes |
|---|---|
| GSE229658 metformin, 388 uM and 776 uM, 48 h | Paired tests across 10 donors; 657 and 3,499 genes, respectively, at BH-adjusted p < 0.05. The 16, 39, and 78 uM contrasts had no genes at that cutoff, so no empty files were created. GEO sample titles say myoblasts, while the series summary describes myotubes. |
| GSE40789 GW501516 (Cardarine), 10 nM, 96 h | Donor-blocked limma analysis of the four paired samples; 6 genes at BH-adjusted p < 0.05, including increased CD36. This is Cardarine, not a separate "Carterin" accession. Reproduction requires R and Bioconductor `limma`. |
| GSE107998 simvastatin and rosuvastatin, 5 uM | GEO's DESeq tables; exported genes with adjusted p < 0.05 (1,805 and 68 genes). CD36 itself did not pass that cutoff (simvastatin adjusted p = 0.059; rosuvastatin adjusted p = 1). |

`results/geo_myotube_signatures/manifest.csv` records contrasts, gene counts, CD36 statistics, methods, and GEO source URLs. The downloaded processed GEO tables and GPL6244 annotation are in `data/geo_myotube/`.

GSE44051 is an electrical-pulse exercise model, not a drug perturbation, and was not added to this drug-signature directory. Its GEO data and prior ranked analysis remain under `data/GSE44051_invitro_stimulation/` and `results/celltype_split/`.

To regenerate the files, run `python notebooks/06_geo_myotube_signatures.py` in the analysis environment with its Python dependencies and R/limma available. The analysis uses only unambiguous Ensembl-to-approved-symbol mappings for GSE229658 and single-symbol platform mappings for GSE40789.
