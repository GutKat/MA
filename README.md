# Master Thesis: RNA Design for Exoribonuclease Resistance (xrRNA) Based on Flavivirus Family

This repository contains the files and code associated with my Master’s thesis, where the goal is to design exoribonuclease-resistant RNA (xrRNA) based on the flavivirus family. The project focuses on two specific flavivirus subfamilies: **MBFV (Mosquito-borne flaviviruses)** and **TBFV (Tick-borne flaviviruses)**.

## Project Structure

### 1. xrRNA Design
The RNA design approaches are besed on two different flavivirus subfamilies: **MBFV** and **TBFV**.

#### Sequence Analysis
- **Path:** `xrRNA_design/analysis_MBFV` and `xrRNA_design/analysis_TBFV`
- **Description:** These folders contain data and Jupyter notebooks used to analyze the sequences of the respective families. It cointains analyse of the structural characteristics of the xrRNA structure. These analyses form the scaffold for designing approach.

#### RNA Design
- **Path:** `xrRNA_design/design_MBFV` and `xrRNA_design/design_TBFV`
- **Important files:**
  - `design.py`: The main script used to design the RNA sequences. This is the most critical file of the project. To execute the design process, run:
    ```bash
    python design.py
    ```
    This starts the design process and outputs a seqeuence.
    This will generate the designed RNA sequence as output.
  - `utils.py`: Contains important functions for the design process, like `mc_optimize()`.
  - `ir_utils.py`: Contains functions for setting up the Infrared model for the RNA process.

  
- **Additional files in MBFV only:** 
  In the `xrRNA_design/design_MBFV` folder, you will also find other files used to generate multiple sequences with (`create_sample.py`) and without optimization (`create_sample_wo_opt.py`). These scripts need an output file or folder, to store the created sequences within the specified file/folder. The output path can be specified by using the flag -o or --output. 
  ```bash
    python create_sample.py -o path_to_folder
  ```
  where the -o defines the path to the output folder and the sequences are stored as 'designs_before_after_opt_i.csv', where i is replaced by a number.
  
    ```bash
      python create_sample_wo_sample.py -o csv_file.csv
    ```
  where the -o defines the path to the output csv file, where the sequences are stored.
  There is also a `leader_design` folder for creating a leader sequence needed for experimental testing.

### 2. Covariance Model
- **Path:** `cov_model/MBFV` and `cov_model/TBFV`
- **Description:** This section contains files for building covariance models based on the MBFV and TBFV sequences.

### 3. Thesis Analysis Notebooks
- **Path:** `thesis`
- **Description:** This folder contains key Jupyter notebooks used for different types of analysis in the thesis.
  - `diversity_analysis.ipynb`: Analyzes the diversity of the RNA sequences of the difference model approaches (Infrared vs. Covariance model)
  - `sequence_analysis.ipynb`: Analyzes the designed RNA sequences.
  - This folder also contains important images used in the thesis.

### 4. Scripts
- **Path:** `scripts`
- **Description:** This folder contains various Python scripts used to automate processes during the development of the thesis.

## Usage Instructions

To run the design process for the RNA sequences, navigate to the appropriate design folder (`xrRNA_design/design_MBFV` or `xrRNA_design/design_TBFV`) and execute the following command:

```bash
python design.py
