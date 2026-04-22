# Neuro-Transcriptomics: Automated Estrogenic Response Pipeline

## Project Overview
This repository contains a specialized bioinformatics pipeline designed to interrogate the human transcriptomic response to 17β-estradiol (E2). By utilizing single-cell RNA-sequencing (scRNA-seq) data from human iPSC-derived neurons, this project models the molecular shifts within the **estrogen-dopamine axis.**

The goal of this investigation is to map the up-regulation of genes involved in cellular growth (IGFBP4/5), metabolic expansion (MYC), and signal transduction, providing a molecular framework for understanding hormonal fluctuations in neuro-biotechnology.

## Technical Features
- **Automated DGE Analysis:** Utilizes Python (Pandas/SciPy) to calculate Log2 Fold Change and P-values for 19,098 genes.
- **API Integration:** Built-in connection to the **MyGene.info API** for real-time automated functional annotation of top-responding genes.
- **Visual Analytics:** Implementation of publication-quality **Volcano Plots** and **Clustered Heatmaps** using Seaborn and Matplotlib.
- **Reproducibility:** Optimized for high-resolution single-cell data (GSE154873).

## Pipeline Architecture
1. **Calibration:** Handshake protocol between expression matrices and clinical metadata.
2. **Analysis:** Statistical thresholding ($p < 0.05$, $|Log2FC| > 1.0$).
3. **Interrogation:** Automated database queries for biological summaries.
4. **Visualization:** Mapping of the "Transcriptomic Explosion" and cellular fingerprints.

## Repository Structure
- `neuro_pipeline.py`: The master execution script.
- `requirements.txt`: Environment configuration.
- `results/`: Contains the generated Volcano Plots and Heatmaps.

## Academic Context
Developed as a Minor Project in Biotechnology (1st Year). This project bridges the gap between computational automation and molecular endocrinology.
