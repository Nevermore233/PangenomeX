# PGcnv: Pangenome-based Incremental Learning Framework for Shallow Copy Number Variation Detection

## Overview
PGcnv is a framework for detecting shallow copy number variations (CNVs) using a pangenome-based incremental learning approach. This repository provides tools for data preprocessing, CNV pre-detection with a CUSUM control chart (ZIP-Caller), and subsequent generation of bubbles, phylogenetic trees, and CNV relationship networks.

## Abstract

## Installation
Uncompress the installation package using the following commands:

```bash
cd /my/install/dir/
unzip /path/to/PGcnv.zip
```

**Requirements**

You'll need to install the following Python packages in order to run the code:

```bash
Python 3.4 or higher
pysam
json
tqdm
argparse
pandas
torch
networkx
numpy
scipy
```

Before starting the project, you need to configure the parameter file. For detailed instructions, refer to **data/parameter_cfg.config**.


## Step 1: Data Preprocessing
Both ZIP-Caller and PGcnv require [data preprocessing] before analysis. The preprocessing consists of two steps: 1.data normalization and 2.baseline setting.

**Usage**
```bash
python3 data_processing.py [-config CONFIG_FILE]

commands:
-config [str]: Path to the configuration file.
```

Example:
```bash
python3 data_processing.py -config data/parameter.config
```
Note: Setting the baseline requires at least 50 normal samples; otherwise, a warning will be issued.

## Step 2: CNV Pre-detection with CUSUM Control Chart (ZIP-Caller)
ZIP-CALLER employs a CUSUM control chart-based model for CNV pre-detection.

**Usage**
```bash
python3 zip-caller.py [-config CONFIG_FILE] [-o OUTPUT_TSV_FILE] [-w SLIDING_WINDOW_SIZE] [-k REFERENCE_VALUE]

commands:
-config [str]: Path to the configuration file.
-o [str]: Path to the output file.
-w [int]: Sliding window size (default: 1000).
-k [float]: Reference value for the allowed degree of deviation (default: 0.3).
```

Example:
```bash
python3 zip_caller.py -config parameter.config -o data/zipcall_output.tsv -w 1000 -k 0.3
```

## Step 3: Generate Bubbles, Phylogenetic Tree, and CNV Relationship Network
Use gen_bubbles.py to generate bubbles for all samples. The bubbles are labeled as 0, 1, or unknown, corresponding to population CNVs, pathogenic CNVs, and pending CNVs, respectively.

**Usage**
```bash
python3 gen_bubbles.py [-config CONFIG_FILE] [-cnv CNV_CSV_FILE] [-o OUTPUT_NPZ_FILE]

commands:
-config [str]: Path to the configuration file.
-cnv [str]: Path to the CNV list in CSV format. 
-o [str]: Path to the output .npz file (default: bub_results.npz).
```
The CSV list file  should be converted from the VCF file of CNV samples in the training dataset. See [example/cnv_list.csv] for details. The CSV should include the following columns: 
[file_name, chr_name, start_pos, end_pos, logr]

Example:
```bash
python3 gen_bubbles.py -config data/parameter.config -cnv data/cnv_list.csv -o data/bub_results.npz
```






