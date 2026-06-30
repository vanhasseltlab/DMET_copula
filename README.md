# Joint distribution modelling of hepatic drug metabolizing enzymes and transporters expression for PBPK-driven virtual patient simulation

These scripts were used to study the joint distribution between the gene expression of hepatic drug metabolizing enzymes and transporters (DMET) to support realistic PBPK predictions which utilize such DMET distributions. 

The project includes the following script files:

-   **01 GTEX processing.R**, processes the GTEx database for downstream analyses.
-   **02 CYP model.R**, develops and evaluates the joint distribution models for the CYP gene family, illustrating the modelling workflow used across all datasets (UGT, SLC, ABC).
-   **03 PBPK-selegiline.R**, investigates the impact of CYP gene joint distribution on PBPK predictions, using selegiline as a case study.
-   **04 decision tree.R**, implements the decision tree used to select appropriate joint distribution models for any dataset of interest.
