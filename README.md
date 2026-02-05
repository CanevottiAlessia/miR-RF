# miR-RF

**miR-RF** is a machine learning framework, composed of 3 modules (RNAfold, miR_application.py and miR_classes.py), for the structural evaluation and classification of human pre-microRNAs, based on features derived from RNA secondary structure predictions.
The software implements the methodology described in the accompanying manuscript (*"An operational workflow for the systematic annotation of human miRNAs"*), providing a reproducible pipeline for pre-miRNA validation and robustness assessment.

---

## Installation

To run miR-RF, you need the following dependencies:
- [Conda](https://docs.conda.io/)
- [Vienna RNAfold](https://anaconda.org/bioconda/viennarna/files?version=cf201901) - required for RNA secondary structure prediction

First, create and configure the Conda environment using the provided configuration file:

```bash
conda env create -f miR_configuration_file.yml
conda activate miR_RF
```

Then, install the additional required Python packages within the miR_RF environment:

```bash
conda install scipy
conda install statsmodels
```

## Usage

After setting up the environment, start from an input FASTA file (<FASTA_file>) and predict RNA secondary structures using RNAfold:

```bash
RNAfold -p -d2 --noLP --noDP --noPS --jobs=<n of threads> <FASTA_file> > <output_RNAfold>
```

Once the RNA secondary structure predictions are generated (<output_RNAfold>), use this file as input for miR_application.py to obtain miR-RF predictions:

```bash
python3 miR_application.py <output_RNAfold> FASTA_file
```
**Output miR_application.py**: "valid" = 2, "non-valid" = 1

Optionally, if you also want to compute the structural stability classes (R, D, I, S), run miR_classes.py using the RNAfold output, the input FASTA file, and the miR-RF predictions:

```bash
python3 miR_classes.py <RNAfold_file> <FASTA_file> <miR-RF_output> <output_file_miR-RF>
```
**Output miR_classes.py**: structural stability class = "status" (R, D, I, S)

---

## Example

FASTA (or multi-FASTA) file as first input: 

```plaintext
>hsa-let-7a-1
UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA

```

Sample input file structure (RNAfold output):

```plaintext
>hsa-let-7a-1
UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA
(((((.(((((((((((((((((((((.....(((...((((....)))).))))))))))))))))))))))))))))) (-34.20)
{((((.(((((((((((((((((((((.....(((...((({....}))).))))))))))))))))))))))))))))} [-35.18]
(((((.(((((((((((((((((((((.....(((...((((....)))).))))))))))))))))))))))))))))) {-34.20 d=3.42}
 frequency of mfe structure in ensemble 0.203686; ensemble diversity 5.63
...
```

miR-RF output

```plaintext
"miRNA name"        "prediction"
">hsa-let-7a-1"     "2"
```

miR-RF_classes output

```plaintext
"miRNA name"        "status"
">hsa-let-7a-1"     "R"
```

---

## Related resource

An interactive web application for exploration and filtering of the annotations generated in the manuscript is available at: https://app-mir-rf-vfd7s8nncj3mx6anbaaxrh.streamlit.app/

---

## Software availability

miR-RF is freely available as an open-source software package at 
https://github.com/CanevottiAlessia/miR-RF.

The interactive miR-RF Annotation Browser is accessible at  
https://app-mir-rf-vfd7s8nncj3mx6anbaaxrh.streamlit.app/  
(source code: https://github.com/CanevottiAlessia/Streamlit-miR-RF).

---

## Citation

If you use miR-RF, please cite: "..."
