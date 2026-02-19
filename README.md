# miR-RF
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
[![Conda](https://img.shields.io/badge/environment-conda-green)](https://docs.conda.io/)


This is the main GitHub repository for **miR-RF**, a machine learning tool and computational workflow for the structural evaluation and classification of human pre-microRNAs. The workflow is composed of 3 main modules: *RNAfold*, *miR_application.py* and *miR_classes.py* and implements key concepts and analytical methods described in the accompanying manuscript (*"An operational workflow for the systematic annotation of human miRNAs"*). 

---

## Installation

miR-RF relies on the Conda virtual environment manager for dependency management and installation.
Please make sure that **Conda** is properly installed on your system before running miR-RF.
- [Conda documentation](https://docs.conda.io/) 

To install miR-RF, first create and configure the Conda environment using the configuration file included in this repository.

```bash
conda env create -f miR_configuration_file.yml
conda activate miR_RF
```

## Usage

**miR-RF** requires as input a FASTA file containing the sequences to be analyzed.

Starting from the input sequences, the workflow proceeds as follows:

1. **RNAfold** is used to predict RNA secondary structures.
2. **miR_application.py** classifies each sequence as a **valid (2)** or **non-valid (1)** pre-miRNA.
3. **miR_classes.py** derives structural stability classes based on features extracted from the predicted secondary structures.

After setting up the Conda environment, start from an input FASTA file (`<FASTA_file>`) and predict RNA secondary structures using **RNAfold**, as follow:

```bash
RNAfold -p -d2 --noLP --noDP --noPS --jobs=<n of threads> <FASTA_file> > <output_RNAfold>
```

For more information on how to use **RNAfold**, including input formats and options, please refer to the official documentation:

- [RNAfold documentation (ViennaRNA Package)](https://viennarna.readthedocs.io/en/latest/tutorial/RNAfold.html)


Once the RNA secondary structure predictions are generated (`<output_RNAfold>`), provide this file as input to **miR_application.py** to obtain **miR-RF** predictions:

```bash
python3 miR_application.py <output_RNAfold> <miR-RF_output>
```
**Output miR_application.py**: the output consists in a tab delineated file with two columns, reporting each sequence in the input and its classification: "valid" = 2 or "non-valid" = 1.

Optionally, if you also want to compute the structural stability classes (R, D, I, S), run **miR_classes.py** using the RNAfold output, the input FASTA file, and miR-RF predictions:

```bash
python3 miR_classes.py <RNAfold_file> <FASTA_file> <miR-RF_output> <output_file_miR-RF>
```

**Output miR_classes.py**: structural stability class = "status" (R, D, I, S). Also in this case the output is provided in the form of a tab delineated table, with one sequence per line and its structural stability class.

--- 

> ⚠️ **Important notes**
>
> **miR-RF does not support FASTA headers containing tab characters**  
> Please ensure that sequence headers do **not** include tabs (`\t`), as this may lead to incorrect parsing or unexpected behavior.  
> We recommend using only whitespace-free identifiers or safe delimiters such as underscores (`_`), pipes (`|`), or semicolons (`;`) in FASTA headers.
>
> **Output and tmp files**  
> Output files from both `miR_application.py` and `miR_classes.py` are written to the directory from which the scripts are launched.  
> Temporary folders and files are created during execution. Temporary files are automatically deleted at the end of the run, while the folders remain empty.


---

## Example

Below you can find a complete example of the workflow. 
An example FASTA file (fasta_example.fa) is provided with the repository (dir: input):

```plaintext
>hsa-mir-3665
GCGGGCGGCGGCGGCGGCAGCAGCAGCAGGUGCGGGGCGGCGGCCGCGCUGGCCGCUCGACUCCGCAGCUGCUCGUUCUGCUUCUCCAGCUUGCGCACCAGCUCC
>hsa-mir-3666
AGUAAGGUCCGUCAGUUGUAAUGAGACCCAGUGCAAGUGUAGAUGCCGACUCCGUGGCAGAGUUCAGCGUUUCACACUGCCUGGUCUCUGUCACUCUAUUGAAUUAGAUUG
>hsa-mir-3667
UGAGGAUGAAAGACCCAUUGAGGAGAAGGUUCUGCUGGCUGAGAACCUUCCUCUCCAUGGGUCUUUCAUCCUCA
>hsa-mir-3668
AUAUAUGAAAUGUAGAGAUUGAUCAAAAUAGUUUCUAUCAAAAUAGUUUUGAUCAAUCUCUGCAAUUUUAUAUAU
>hsa-mir-3670-1
UCUAGACUGGUAUAGCUGCUUUUGGAGCCUCACCUGCUGAGAGCUCACAGCUGUCCUUCUCUAGA
```
Users may replace this file with any FASTA or multi-FASTA file of interest.

**RNA secondary structure prediction**

RNA secondary structures are predicted using RNAfold.
The name of the output file can set by the user (here RNAfold_file.txt):

```bash
RNAfold -p -d2 --noLP --noDP --noPS --jobs=<n of threads> fasta_example.fa > RNAfold_example.txt
```

Example output (RNAfold_example.txt):

```plaintext
>hsa-mir-3665
GCGGGCGGCGGCGGCGGCAGCAGCAGCAGGUGCGGGGCGGCGGCCGCGCUGGCCGCUCGACUCCGCAGCUGCUCGUUCUGCUUCUCCAGCUUGCGCACCAGCUCC
(((((((..(((((.(((((.((((((((.(((((((((((((((.....)))))).)).))))))).))))).))))))))...)).)))..))).)).))... (-64.40)
(((((((..(((((.(((((.((((((((.(((((((((((((((,...,)))))).)).))))))).))))).))))))))...)).)))..))).)).))... [-65.07]
(((((((..(((((.(((((.((((((((.(((((((((((((((.....)))))).)).))))))).))))).))))))))...)).)))..))).)).))... {-64.40 d=3.20}
 frequency of mfe structure in ensemble 0.334897; ensemble diversity 5.75
>hsa-mir-3666
AGUAAGGUCCGUCAGUUGUAAUGAGACCCAGUGCAAGUGUAGAUGCCGACUCCGUGGCAGAGUUCAGCGUUUCACACUGCCUGGUCUCUGUCACUCUAUUGAAUUAGAUUG
.....((((..(((((.((...((((((.((.(((.((((((((((.(((((.......)))))..)))))).)))))))))))))))....))...)))))....)))). (-32.30)
.....,(((..(((((.({...(((({{{((.(({{((((((((((.(((((.......))))}}.)))))).)))))))))))))))..,.)}...)))))....))),. [-35.43]
......(((..(((((.((...(((((.(((.((..((((((((((.(((((.......)))).).)))))).)))).))))))))))....))...)))))....))).. {-29.10 d=11.82}
 frequency of mfe structure in ensemble 0.0062276; ensemble diversity 17.50
>hsa-mir-3667
UGAGGAUGAAAGACCCAUUGAGGAGAAGGUUCUGCUGGCUGAGAACCUUCCUCUCCAUGGGUCUUUCAUCCUCA
((((((((((((((((((.((((.(((((((((........))))))))).)))).)))))))))))))))))) (-51.60)
{(((((((((((((((((.((((.(((((((((........))))))))).)))).)))))))))))))))))} [-52.19]
((((((((((((((((((.((((.(((((((((........))))))))).)))).)))))))))))))))))) {-51.60 d=1.05}
 frequency of mfe structure in ensemble 0.38599; ensemble diversity 1.55
>hsa-mir-3668
AUAUAUGAAAUGUAGAGAUUGAUCAAAAUAGUUUCUAUCAAAAUAGUUUUGAUCAAUCUCUGCAAUUUUAUAUAU
(((((((((((((((((((((((((((((.((((......)))).)))))))))))))))))).))))))))))) (-33.70)
{((((((((((((((((((((((((((((.((((......)))).))))))))))))))))))},)))))))))} [-34.64]
(((((((((((((((((((((((((((((.((((......)))).))))))))))))))))))).)))))))))) {-33.70 d=1.93}
 frequency of mfe structure in ensemble 0.216196; ensemble diversity 2.55
>hsa-mir-3670-1
UCUAGACUGGUAUAGCUGCUUUUGGAGCCUCACCUGCUGAGAGCUCACAGCUGUCCUUCUCUAGA
((((((..((.(((((((......((((((((.....)))).)))).)))))))))...)))))) (-24.20)
{(((((..((.(((((((....,.((((((((.....)))).)))),)))))))))...)))))} [-25.20]
((((((..((.(((((((......((((((((.....)))).)))).)))))))))...)))))) {-24.20 d=2.82}
 frequency of mfe structure in ensemble 0.197459; ensemble diversity 4.22
```


**miR-RF prediction**

RNAfold output is processed by miR_application.py to predict if -based on sequence and secondary structure- the input sequence qualify as "valid" candidate pre-miRNAs.
The name of the output file (here "output_miR_application.txt") can be set by the user at runtime:

```bash
python3 miR_application.py RNAfold_example.txt output_miR_application.txt
```

Example output (output_miR_application.txt):

```plaintext
"miRNA name"    "prediction"
">hsa-mir-3665" "1"
">hsa-mir-3666" "1"
">hsa-mir-3667" "2"
">hsa-mir-3668" "2"
">hsa-mir-3670-1" "2"
```


**Structural stability class assignment (optional)**

Finally, miR_classes.py can assign structural stability classes (R, D, I, S).
Also in this case, users are free to set the name of the output file.

```bash
python3 miR_classes.py RNAfold_example.txt fasta_example.fa output_miR_application.txt output_miR-RF_classes.txt
```

Example output (output_miR-RF_classes.txt):

```plaintext
mirna_name      new_status
>hsa-mir-3665 S
>hsa-mir-3666 I
>hsa-mir-3667 R
>hsa-mir-3668 R
>hsa-mir-3670-1 D
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

## License

This work is licensed under a Creative Commons Attribution 4.0 International License (CC BY 4.0).

You are free to:
- Share — copy and redistribute the material in any medium or format
- Adapt — remix, transform, and build upon the material for any purpose, even commercially

Under the following terms:
- Attribution — appropriate credit must be given to the original authors and the accompanying manuscript.

For details, see: https://creativecommons.org/licenses/by/4.0/

---

## Citation

If you use miR-RF, please cite:

Canevotti A. et al. *An operational workflow for the systematic annotation of human miRNAs*.  
Manuscript under peer review.
