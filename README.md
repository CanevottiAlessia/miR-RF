# miR-RF

**miR-RF** is a machine learning framework, composed of 3 main modules (RNAfold, miR_application.py and miR_classes.py), for the structural evaluation and classification of human pre-microRNAs, based on features derived from RNA secondary structure predictions.
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
python3 miR_application.py <output_RNAfold> <miR-RF_output>
```
**Output miR_application.py**: "valid" = 2, "non-valid" = 1

Optionally, if you also want to compute the structural stability classes (R, D, I, S), run miR_classes.py using the RNAfold output, the input FASTA file, and the miR-RF predictions:

```bash
python3 miR_classes.py <RNAfold_file> <FASTA_file> <miR-RF_output> <output_file_miR-RF>
```

**Output miR_classes.py**: structural stability class = "status" (R, D, I, S)

---

## Example

Below we report a complete example workflow using a FASTA (or multi-FASTA) input file.

Input FASTA file

An example FASTA file (example_FASTA_file.fa) is provided with the repository:

```plaintext
>hsa-let-7a-3_MI0000062
GGGUGAGGUAGUAGGUUGUAUAGUUUGGGGCUCUGCCCUGCUAUGGGAUAACUAUACAAUCUACUGUCUUUCCU
>hsa-let-7b_MI0000063
CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUGUUGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG
>hsa-let-7c_MI0000064
GCAUCCGGGUUGAGGUAGUAGGUUGUAUGGUUUAGAGUUACACCCUGGGAGUUAACUGUACAACCUUCUAGCUUUCCUUGGAGC
```
Users may replace this file with any custom FASTA or multi-FASTA file of interest.

RNA secondary structure prediction

RNA secondary structures are predicted using RNAfold.
The output file name can be freely chosen by the user (here named RNAfold_file.txt):

```bash
RNAfold -p -d2 --noLP --noDP --noPS --jobs=<n of threads> example_FASTA_file.fa > RNAfold_file.txt
```

Example output (RNAfold_file.txt):

```plaintext
>hsa-let-7a-3_MI0000062
GGGUGAGGUAGUAGGUUGUAUAGUUUGGGGCUCUGCCCUGCUAUGGGAUAACUAUACAAUCUACUGUCUUUCCU
(((.(((((((((((((((((((((((((((...)))))).........))))))))))))))))))))).))) (-34.10)
(((.(((((((((((((((((((((((((((...)))))).........))))))))))))))))))))).))} [-34.81]
(((.(((((((((((((((((((((((((((...)))))).........))))))))))))))))))))).))) {-34.10 d=2.00}
 frequency of mfe structure in ensemble 0.314643; ensemble diversity 3.41
>hsa-let-7b_MI0000063
CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUGUUGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG
(((((.(((((((((((((((((((((((.((((((.....))))))...))).....))))))))))))))))))))))))) (-46.70)
(((((.(((((((((((((((((((((((,((((((.....)))))).,.,}}....}))))))))))))))))))))))))) [-48.38]
(((((.(((((((((((((((((((((...((((((.....))))))..........)))))))))))))))))))))))))) {-46.20 d=5.25}
 frequency of mfe structure in ensemble 0.0651573; ensemble diversity 7.60
>hsa-let-7c_MI0000064
GCAUCCGGGUUGAGGUAGUAGGUUGUAUGGUUUAGAGUUACACCCUGGGAGUUAACUGUACAACCUUCUAGCUUUCCUUGGAGC
((.((((((..(((.(((.(((((((((((((..((.((.((...)).)).))))))))))))))).))).)))..)))))))) (-31.40)
((.((((((..(((.(((.(((((((((((((..((.(,.({...}).,).))))))))))))))).))).)))..)))))))) [-33.10]
((.((((((..(((.(((.(((((((((((((..((.(..(.....)..).))))))))))))))).))).)))..)))))))) {-31.20 d=4.73}
 frequency of mfe structure in ensemble 0.0635097; ensemble diversity 7.33
```

miR-RF prediction

The RNAfold output is then used as input for miR_application.py to obtain miR-RF predictions.
The output file name can be chosen arbitrarily (here output_miR_application.txt):

```bash
python3 miR_application.py RNAfold_file.txt output_miR_application.txt
```

Example output (output_miR_application.txt):

```plaintext
"miRNA name"        "prediction"
">hsa-let-7a-3_MI0000062"     "2"
">hsa-let-7b_MI0000063"     "2"    
">hsa-let-7c_MI0000064"     "2"
```

Structural stability class assignment (optional)

Finally, miR_classes.py can be run to assign structural stability classes (R, D, I, S).
All output file names can be freely specified by the user:

```bash
python3 miR_classes.py example_RNAfold_file.txt example_FASTA_file.fa output_miR_application.txt output_miR-RF_classes.txt
```

Example output (output_miR-RF_classes.txt):

```plaintext
"miRNA name"        "status
">hsa-let-7a-3_MI0000062"     "R"
">hsa-let-7b_MI0000063"     "R"    
">hsa-let-7c_MI0000064"     "R"
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
