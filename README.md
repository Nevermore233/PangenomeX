# PGcnv: Pangenome-based Incremental Learning Framework for Shallow Copy Number Variation Detection

## Abstract
Shallow whole-genome sequencing (sWGS) for inferring copy number variations (CNVs) in the human genome has rapidly become a routine practice at many research centers. Although several tools are available for inferring CNVs from sWGS data, the precision of these methods is often limited due to noise from low coverage and the impact of copy number polymorphisms (CNPs), which makes them insufficient for large-scale variant detection in clinical practice.In this study, we developed PGcnv, a tool that includes a CNV caller (ZIP-Caller) based on a zero-inflated Poisson distribution and a corrector based on a pangenome graph with bubbles. First, we develop a detection method specifically designed to account for the depth distribution characteristics of sWGS data, ensuring high sensitivity. Next, we incorporate a CNV relationship network that integrates and learns CNPs from the population. Using a graph convolutional network model, PGcnv effectively differentiates between noise and true copy number variations. Additionally, it supports incremental learning, allowing the model to improve as more data is incorporated. We tested PGcnv on simulated datasets and 561 real-world samples from four major regions of China, and compared its performance with six other tools. The results demonstrate that PGcnv outperforms current popular tools. Overall, offering a novel perspective for population-based CNV detection. 


## Installation
Uncompress the installation package using the following commands:

```bash
cd /my/install/dir/
unzip /path/to/PGcnv.zip
```

**Requirements**

You'll need to install the following Python packages in order to run the code:

```bash
Python 3.11.9 
pysam 0.22.1
json 1.9.6
tqdm 4.67.1
pandas 2.2.3
torch 2.6.0
torchaudio 2.6.0             
torchvision 0.21.0 
networkx 3.4.2
numpy 2.2.3
scipy 1.15.2
```

Before starting the project, you need to configure the parameter file. For detailed instructions, refer to **my.config**. 


## Step 1: Data Preprocessing
Both ZIP-Caller and PGcnv require [data preprocessing] before analysis. The preprocessing consists of two steps: 1.data normalization and 2.baseline setting.

**Usage**
```bash
python3 data_processing.py [-config CONFIG]

commands:
-config [str]: Path to the configuration file.
```

Example:
```bash
python3 data_processing.py -config my.config
```
This step will also output ' rj_means_and_n.json' in the current folder for incremental updates.

Note: Setting the baseline requires at least 50 normal samples; otherwise, a warning will be issued.

## Step 2: CNV Pre-detection with CUSUM Control Chart (ZIP-Caller)
ZIP-Caller employs a CUSUM control chart-based model for CNV pre-detection.

**Usage**
```bash
python3 zip-caller.py [-config CONFIG] [-o OUTPUT] [-w SLIDING_WINDOW_SIZE] [-k REFERENCE_VALUE]

commands:
-config [str]: Path to the configuration file.
-o [str]: Path to the output file.
-w [int]: Sliding window size (default: 3000).
-k [float]: Reference value for the allowed degree of deviation (default: 0.3).
```

Example:
```bash
python3 zip_caller.py -config my.config -o data/zipcall-output -w 3000 -k 0.3
```

## Step 3: Generate Bubbles, Phylogenetic Tree, and CNV Relationship Network
Use gen_bubbles.py to generate bubbles for all samples. The bubbles are labeled as 0, 1, or unknown, corresponding to population CNVs, pathogenic CNVs, and pending CNVs, respectively.

**Usage**
```bash
python3 gen_bubbles.py [-config CONFIG] [-cnv CSV] [-o NPZ]

commands:
-config [str]: Path to the configuration file.
-clist [str]: Path to the CNV list in CSV format. 
-o [str]: Path to the output .npz file (default: bub_results.npz).
```
The CSV list file  should be converted from the VCF file of CNV samples in the training dataset. See [input_csv/cnv_list.csv] for details. The CSV should include the following columns: 
[file_name, chr_name, start_pos, end_pos, logr]

Example:
```bash
python3 gen_bubbles.py -config my.config -clist input_csv/cnv_list.csv -o bub_results.npz
```
Then, use tree2graph.py to convert the .nwk file to a .gpickle file.

**Usage**
```bash
python3 tree2graph.py [-nwk NWK] [-npz NPZ] [-cnv CNV] [-o GPICKLE]

commands:
-nwk [str]: Path to the .nwk file.
-npz [str]:Path to the bubble_result.npz file(the output of gen_bubbles.py).
-cnv [str] Path to the .cnv file(the output of zip_caller).
-o [str]: Path to the output(.gpickle) file, example: graph.gpickle
-update [store_true]: When this parameter appears in the command line, it indicates that the tree2graph.py will use the updated files as input (see Incremental Update Module).

```
Example:
```bash
python3 tree2graph.py -nwk mynwk.nwk -npz bub_results.npz -cnv data/zipcall-output/zipcaller_res_2025-04-01_16-25-47.cnv -o graph.gpickle
```
This step will also output 'df_edge.csv ' in the current folder for incremental updates. 

## Step 4: CNV calling based on graph convolutional network
Use PGcnv.py to obtain the final detection results.

**Usage**
```bash
python3 pgcnv.py [-config CONFIG] [-k GPICKLE] [-o PATH]

commands:
-config [str]: Path to the '.config' file.
-k [str]: Path to the '.gpickle' file
-o [str]: Path to the output file, example: data/output

```
Example:
```bash
python3 pgcnv.py -config my.config -k graph.gpickle -o data/pgcnv_output
```

## Incremental Update Module
This module achieves incremental learning by updating the normalization file and the CNV relationship network.Use 'incremental_update.py' to obtain the updated files.

```bash
python3 Incremental_update.py [-config CONFIG] [-i CSV] [-rj JSON] [-bub NPZ] [-edge CSV]

commands:
-config [str]: Path to the configuration file
-i [str]: Path to the 'new_df.csv' file
-rj [str]: Path to the 'rj_means_and_n.json' file
-bub [str]: Path to the 'bub_results.npz' file
-edge [str]: Path to the 'df_edge.csv' file
```
Where rj_means_and_n.json is the file generated by data_processing.py, and df_edge.csv is an intermediate result from tree2graph.py. These two files are stored in the current directory by default.
File bub_results.npz is the output file from gen_bubbles.py.

Example:
```bash
python3 incremental_update.py -config my.config -rj rj_means_and_n.json -bub data/bub_results.npz -edge df_edge.csv
```

The incremental_update.py file will output the updated normalization files for the new samples in the 'data/nor/' directory under the current folder. Additionally, it will generate 'updated_rj_means_and_n.json', 'new_bub_results.npz', and 'new_df_edge.csv' in the current folder. These three files can be used for the next incremental update.

After obtaining these files, you can use the -update flag in tree2graph to generate the new [.gpickle] file.

Example:
```bash
python3 tree2graph.py -nwk mynwk.nwk -npz bub_results.npz -cnv data/zipcall-output/zipcaller_res_2025-04-01_16-25-47.cnv -o graph.gpickle -update
```

At this point, the program will use the updated files to proceed with this step and then proceed with the steps outlined earlier.
