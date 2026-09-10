# Source data

This release contains two complementary source-data workbooks.

## Base manuscript source data

`HeptadPair_Source_Data_v8.5_BASE_MANUSCRIPT_20260830.xlsx` preserves the pre-Phase81/82 source-data package underlying the established CCNG1, CCmax, orthogonal-recovery, CC0, MP3 and supplementary stress-test results. Its 16 sheets include the GLMY/path-complex readout audit. Three cells in `Feature_Accounting` are simple cumulative formulas; they are retained from the frozen workbook.

Sheets:

- `Feature_Accounting`
- `Dataset_Provenance`
- `Feature_Contract`
- `CCNG1_Main_Results`
- `CCmax_ZeroShot`
- `CCmax_Locked`
- `Orthogonal_Recovery`
- `CC0_Results`
- `MP3_Results`
- `Supplementary_Potapov`
- `Supplementary_SYNZIP`
- `Supplementary_Structure_Audit`
- `GLMY5_Readout_Audit`
- `GLMY5_Paired_Deltas`
- `GLMY5_Selection_Freq`
- `GLMY5_FiniteState`

## v8.9 Phase81/Phase82 revision source data

`HeptadPair_Source_Data_v8.9_GITHUB_20260903.xlsx` is the public-safe v8.9 revision workbook. Local filesystem prefixes were replaced with repository-relative placeholders; numerical values and sheet organization are unchanged. It contains the full-factorial E/O/L analysis and the modeled structural-case audit used for revised Figures 2 and 4 and Supplementary Figure S6.

Sheets:

- `README`
- `Fig2_Frozen_CCNG1`
- `Fig2_Factorial_Splits`
- `Fig2_Factorial_Summary`
- `Fig2_LOO_Splits`
- `Fig2_LOO_Summary`
- `Fig2_Shapley_Summary`
- `Fig4_Retrieval_Metrics`
- `Fig4_Set_Recovery`
- `Phase82_Case_Selection`
- `Phase82_Offtarget_Candidates`
- `Fig4_Structure_Cases`
- `Phase82_Model_Geometry`
- `Phase82_Paired_Contrasts`
- `FigS6_4H1802_Parallel`

The CSV exports in `csv/` correspond to the Phase81/Phase82 revision tables used by the portable figure-reproduction script.
