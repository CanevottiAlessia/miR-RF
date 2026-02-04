# miR-RF

**miR-RF** is a machine learning framework for the structural evaluation and classification of human pre-microRNAs, based on features derived from RNA secondary structure predictions.
The software implements the methodology described in the accompanying manuscript (*[titolo paper]*), providing a reproducible pipeline for pre-miRNA validation and robustness assessment.

---

## Summary

miR-RF integrates RNA secondary structure features with a Random Forest classifier to:
- assess pre-miRNA structural viability;
- evaluate sensitivity to single-nucleotide variation (SNV);
- classify pre-miRNAs into biologically interpretable robustness classes.

---

## Components

The repository includes two command-line tools:

| Tool | Function |
|------|----------|
| **miR-RF** | Binary classification of pre-miRNA structural viability |
| **miR-RF_classes** | Classification of pre-miRNAs into robustness classes (R, D, I, S) based on *in silico* SNV analysis |

---

## miR-RF

**miR-RF** predicts whether a pre-miRNA is structurally compatible with canonical miRNA biogenesis.

- **Input:** RNAfold output generated with [ViennaRNA RNAfold](https://www.tbi.univie.ac.at/RNA/)
- **Output:** Viability label ("valid" = 2, "non-valid" = 1)

Predictions are based exclusively on RNA secondary structure features and do not rely on expression or evolutionary conservation data.

---

## miR-RF_classes

**miR-RF_classes** evaluates the robustness of pre-miRNAs to sequence variation by simulating all possible single-nucleotide variants and comparing predictions to the reference sequence.

Each pre-miRNA is assigned to one of four classes:

- R (Resilient): pre-miRNAs initially classified as "valid" not showing a significant enrichment of LoF variants;
- D (Dispensable): pre-miRNAs initially classified as "valid" showing a significant enrichment of LoF variants;
- I (Inducible): pre-miRNAs initially classified as "non-valid" showing a significant enrichment of GoF variants;
- S (Spurious): pre-miRNAs initially classified as "non-valid" not showing a significant enrichment of GoF variants.

---

## Dependencies

- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) (RNAfold)
- [Conda](https://docs.conda.io/)

---

## Installation

Dependencies are managed via Conda:

```bash
conda env create -f miR_configuration_file.yml
conda activate miR_RF
```

ViennaRNA (RNAfold) is required for structure prediction.

Also install scipy and statsmodels with conda in miR-RF environment (conda install scipy, conda install statsmodels)

## Usage

Run miR-RF:

```bash
python3 miR_application.py <RNAfold_input> <output_file>
```

Run miR-RF_classes:

```bash
python3 miR_classes.py <RNAfold_file> <FASTA_file> <miR-RF_output> <output_file>
```

---

## Related resource

An interactive web application for exploration and filtering of the annotations generated in the manuscript is available at: https://app-mir-rf-vfd7s8nncj3mx6anbaaxrh.streamlit.app/

---

## Software availability

miR-RF is freely available as an open-source software package at 
https://github.com/CanevottiAlessia/miR-RF.

The interactive pre-miRNA Annotation Browser is accessible at  
https://app-mir-rf-vfd7s8nncj3mx6anbaaxrh.streamlit.app/  
(source code: https://github.com/CanevottiAlessia/Streamlit-miR-RF).

---

## Citation

If you use miR-RF, please cite: "..."
