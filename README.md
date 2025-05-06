# Sex Differences in Alzheimer's Disease: Single-cell RNA-seq Analysis

This repository contains code for analyzing sex differences in Alzheimer's Disease (AD) at the cellular level using single-cell RNA sequencing data.

## Overview

We investigate sex-based differences in AD using single-cell RNA-seq data from three human datasets, spanning both sexes across control and disease conditions. Our analysis reveals sex-specific transcriptional signatures and cellular responses to AD pathology.

## Repository Structure

### Dataset Folders
- `GSE147528/`, `GSE160936/`, `GSE174367/` - Raw data from three GEO datasets
- `cell_types/` - Cell type annotation files and markers

### Processing Scripts
- `download_data.py` - Downloads datasets from GEO
- `mapping.py` - Maps gene names to Ensembl IDs for consistency
- `merge_data.py` - Combines data from different sources
- `processing_data.py` - Core data preprocessing functions
- `preprocess.r`, `preprocess.ipynb` - R and Python preprocessing workflows
- `convert_all.R` - Converts data formats for compatibility
- `conversion.py` - Python utilities for data conversion
- `check_corrupted_qsave.r` - Validates R data objects
- `check_empty.py` - Identifies empty datasets
- `check_gene.py` - Verifies gene presence across datasets
- `loop_all_ann.py` - Batch processes annotation files

### Analysis Scripts
- `deseq2.py` - Differential expression analysis using PyDESeq2
- `mrvi_train.py` - Implements Multi-Resolution Variational Inference
- `enrichment.py`, `enrichment_2.py` - Gene set enrichment analysis
- `plot_umap.py` - Cell clustering and dimensionality reduction
- `view_data.py` - Data exploration utilities

### Visualization
- `violin_plots.py` - Creates violin plots for gene expression
- Multiple GSEA plot files (`.png`, `.pdf`) - Visualizations of enrichment results
- `custom_gsea_neurogenesis.png` - Custom GSEA visualization for neurogenesis pathways

### Exploratory Analysis
- `eda.ipynb`, `new_eda.ipynb` - Jupyter notebooks for exploratory data analysis
- `1_DE/` - Differential expression results
- `1_Enrichment/` - Enrichment analysis results
- `_process_again/` - Secondary processing workflows

### Data Files
- `all.csv` - Combined dataset
- `common_genes.csv` - Genes found across all datasets
- `corrupted.csv` - List of corrupted data files
- `empty_datasets.csv` - List of empty datasets
- `genes_list.txt` - Reference gene list

### Environment Configuration
- `requirements.txt` - Python package dependencies
- `ssread_env.yml` - Conda environment specification

## Data

The analysis uses three publicly available datasets from GEO:
- GSE147528
- GSE160936
- GSE174367

## Key Features

- Cell typing and gene name mapping to Ensembl IDs
- Differential expression analysis using PyDESeq2
- Multi-Resolution Variational Inference (MrVI) for differential abundance analysis
- Gene set enrichment analysis (GSEA)
- Visualization of sex-specific gene expression patterns

## Key Findings

- Identification of female-elevated protective factors (ANGPTL4, HSPH1, MT2A, KDM6A)
- Male-specific transcriptional signatures (TTTY14, USP9Y, LINC00278, NLGN4Y)
- Sex-dependent pathway activation in major brain cell types

## Usage

1. Set up the environment:
   ```
   conda create env -n ssread
   conda activate ssread
   pip install -r requirements.txt
   ```

2. Download the datasets using `download_data.py`

3. Process the data:
   - Cell typing with `plot_umap.py`
   - Gene mapping with `mapping.py`
   - Merge datasets with `merge_data.py`

4. Run analyses:
   - Differential expression with `deseq2.py`
   - MrVI training with `mrvi_train.py`
   - Enrichment analysis with `enrichment.py`

5. Visualize results with various plotting scripts

## Authors

- Darya Taratynova
- Arman Bolatov
- Shahad Hardan