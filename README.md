# Master Thesis: RNA Design for Exoribonuclease Resistance (xrRNA) Based on Flavivirus Family

This repository contains the files and code associated with my Master’s thesis, where the goal is to design exoribonuclease-resistant RNA (xrRNA) based on the flavivirus family. The project focuses on two specific flavivirus subfamilies: **MBFV (Mosquito-borne flaviviruses)** and **TBFV (Tick-borne flaviviruses)**.

## Project Structure

### 1. xrRNA Design
The main focus of this project is to design RNA sequences that are resistant to exoribonucleases. The sequences are designed for two different flavivirus subfamilies: **MBFV** and **TBFV**.

#### Sequence Analysis
- **Path:** `xrRNA_design/analysis_MBFV` and `xrRNA_design/analysis_TBFV`
- **Description:** These folders contain data and Jupyter notebooks used to analyze the sequences of the respective families. These analyses form the basis for designing the xrRNA.

#### RNA Design
- **Path:** `xrRNA_design/design_MBFV` and `xrRNA_design/design_TBFV`
- **Important files:**
  - `design.py`: The main script used to design the RNA sequences. This is the most critical file of the project. To execute the design process, run:
    ```bash
    python design.py
    ```
    This will generate the designed RNA sequence as output.
  - `utils.py` and `ir_utils.py`: Helper files used for supporting functionalities in the design process.
  
- **Additional files (MBFV only):** 
  In the `xrRNA_design/design_MBFV` folder, you will also find other files used to generate multiple sequences with and without optimization. There is also a `leader_design` folder for creating a leader sequence needed for experimental testing.

### 2. Covariance Model
- **Path:** `cov_model/MBFV` and `cov_model/TBFV`
- **Description:** This section contains files for building covariance models based on the MBFV and TBFV sequences.

### 3. Thesis Analysis Notebooks
- **Path:** `thesis`
- **Description:** This folder contains key Jupyter notebooks used for different types of analysis in the thesis.
  - `diversity_analysis.ipynb`: Analyzes the diversity of the RNA sequences.
  - `sequence_analysis.ipynb`: Analyzes the designed RNA sequences.
  - This folder also contains important images used in the thesis.

### 4. Scripts
- **Path:** `scripts`
- **Description:** This folder contains various Python scripts used to automate processes during the development of the thesis.

## Usage Instructions

To run the design process for the RNA sequences, navigate to the appropriate design folder (`xrRNA_design/design_MBFV` or `xrRNA_design/design_TBFV`) and execute the following command:

```bash
python design.py
