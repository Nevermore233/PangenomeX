# PGcnv: Pangenome-based Incremental Learning Framework for Shallow Copy Number Variation Detection

## Overview
PGcnv is a framework for detecting shallow copy number variations (CNVs) using a pangenome-based incremental learning approach. This repository provides tools for data preprocessing, CNV pre-detection with a CUSUM control chart (ZIP-Caller), and subsequent generation of bubbles, phylogenetic trees, and CNV relationship networks.

## Installation
Uncompress the installation package using the following commands:

```bash
cd /my/install/dir/
unzip /path/to/PGcnv.zip
```

**Requirements**

You'll need to install the following Python packages in order to run the code:

```bash
- Python 3.4 or higher
- pysam
- json
- tqdm
- argparse
- pandas
- torch
- networkx
- numpy
- scipy
```

Before starting the project, you need to configure the parameter file. For detailed instructions, refer to <div style="display: inline-block; padding: 2px 4px; box-shadow: 2px 2px 5px grey;">data/parameter_cfg.config</div>. 

