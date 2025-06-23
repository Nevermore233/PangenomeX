# PangenomeX: a graph convolutional network-based pangenome framework for unbiased population-scale genomic variation analysis

## Abstract
In population-scale genomic variation studies based on shallow whole genome sequencing (sWGS), pangenomes have become a powerful and effective auxiliary tool. They excel at identifying population-specific single nucleotide polymorphisms (SNPs) and insertions/deletions (indels) while mitigating the limitations imposed by shallow read depth. Extending these advantages to copy number variation (CNV), however, remains challenging because two key issues are still unresolved. First, current pangenome frameworks exhibit pronounced population-representation bias arising from uneven sampling across populations. As a result, the graphs primarily capture variants from majority groups while dampening signals from minority groups. Although this bias has limited impact on detecting SNPs and indels, because their accurate identification can often be achieved without explicit population context when sequencing depth is sufficient, it severely compromises CNV detection accuracy, especially under shallow-coverage conditions. Second, in population-scale genomic variation analyses, common but benign population-specific copy number polymorphisms (CNPs) frequently obscure pathogenic CNVs. Existing pangenome frameworks lack dedicated mechanisms for representing CNPs and CNVs, limiting their ability to distinguish pathogenic CNVs from benign, population-specific CNPs. In this study, we present PangenomeX, a graph-convolutional pangenome framework tailored to low-coverage, population-scale CNV analysis. To address the challenge of representing CNPs, we embed known CNPs as prior knowledge into the pangenome graph and construct a CNV relationship network guided by a phylogenetic tree. We then employ a graph convolutional network to learn the interactions between CNV and CNP nodes. To counter population-representation bias, the GCN restricts aggregation to one- and two-hop neighborhoods each node, preserving essential local population context while preventing signals from majority groups from dominating those of minority populations. Evaluation on simulated cohorts and 561 real samples shows that PangenomeX distinguishes pathogenic CNVs from common population CNPs markedly better than existing methods. Overall, PangenomeX offers a methodological blueprint for large-cohort variant screening and provides a practical path for bringing graph-based genomics into clinical practice.


## Installation
Uncompress the installation package using the following commands:

```bash
cd /my/install/dir/
unzip /path/to/PangenomeX.zip
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
Both ZIPcnv and PangenomeX require [data preprocessing] before analysis. The preprocessing consists of two steps: 1.data normalization and 2.baseline setting.

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
