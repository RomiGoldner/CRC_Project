# Colorectal Cancer (CRC) Project

This project includes TCR-sequencing analysis for 3 different CRC datasets, one bulk datasets introduced here for the first time and two published single-cell datasets.

<p align="center">
<img src="plots/Figure 1 CRC Article.png" width="640" height="480">
</p>

## Installation and Usage

To Clone and Set Environment:
1. Clone the GitHub repository and create its requisite conda environment as follows.<br />
   Make sure you use a recent conda version, e.g. version=4.10 or above

```bash
conda env create -n my_env_name --file=environment.yml

conda activate my_env_name
```

2. Some of the analysis is based on embeddings created by the [scCVC model](https://www.science.org/doi/10.1126/sciadv.adk4670).
   Follow the instructions [here](https://github.com/RomiGoldner/CVC) to install and create the embeddings. 

## Notebooks
The main notebooks used in the project are under the notebooks folder and classification folder. <br />

The [notebooks](https://github.com/RomiGoldner/CRC_Project/tree/main/notebooks) folder contains notebooks that are used to do the individual and over analysis on the different datasets. 
- [notebooks/CRC_data_embedding_analysis.ipynb](https://github.com/RomiGoldner/tree/main/CRC_Project/notebooks/CRC_data_embedding_analysis.ipynb) includes the general analysis for the CRC TCR bulk sequencing dataset. <br />
- [notebooks/CRC_low_high_risk_analysis.ipynb](https://github.com/RomiGoldner/tree/main/CRC_Project/notebooks/CRC_low_high_risk_analysis.ipynb) includes the risk analysis. <br />
- [notebooks/GSE164522_general_analysis.ipynb](https://github.com/RomiGoldner/tree/main/CRC_Project/notebooks/GSE164522_general_analysis.ipynb) includes the general analysis for the first single cell sequencing dataset (GSE164522)[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164522]. <br />
- [notebooks/Zhang_general_analysis.ipynb](https://github.com/RomiGoldner/CRC_Project/tree/main/notebooks/Zhang_general_analysis.ipynb) includes the general analysis for the second single cell sequencing (dataset)[http://crctcell.cancer-pku.cn/]. <br />
- [notebooks/all_3_datasets_comparison.ipynb](https://github.com/RomiGoldner/CRC_Project/tree/main/notebooks/all_3_datasets_comparison.ipynb) is the overlap analysis on all 3 datasets. <br />


The [classification](https://github.com/RomiGoldner/CRC_Project/tree/main/classification) folder contains notebooks that are useful for generating the classification results:. <br />
- [MAIT Classification](https://github.com/RomiGoldner/CRC_Project/tree/main/classification/CRC_MAIT_classification.ipynb). <br />
- [Tissue Classification](https://github.com/RomiGoldner/CVC/blob/main/classification/tissue_classfication_2.ipynb). <br />
- Risk Classification is at the bottom of the [risk analysis notebook](https://github.com/RomiGoldner/tree/main/CRC_Project/notebooks/CRC_low_high_risk_analysis.ipynb). <br />
