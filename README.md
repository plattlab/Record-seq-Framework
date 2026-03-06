# An End To End Computational Framework for Record-seq Transcriptional Recording Data

This repository provides the computational framework and benchmarking scripts required to reproduce the results presented in: **"An End-To-End Computational Framework for 'Record-seq' Transcriptional Recording Data"** (doi placeholer).

## Framework Description

Record-seq enables the recording of transcriptional activity in bacterial populations. The technology utilizes an engineered Cas1-Cas2 acquisition complex from *Fusicatenibacter saccharivorans*, in which Cas1 is fused to a reverse transcriptase to capture intracellular RNA fragments. These fragments are integrated as spacers into CRISPR arrays on a recording plasmid, allowing for the reconstruction of past cellular transcriptional histories. Record-seq has been successfully applied to monitor microbial responses to dietary perturbations and inflammatory states within the mouse gut microbiome.

### Previous Publications
* **Florian Schmidt, Maria Y. Cherepkova, and Randall J. Platt.** *Transcriptional recording by CRISPR spacer acquisition from RNA.* Nature, 562(7727):380–385, 2018. (https://www.nature.com/articles/s41586-018-0569-1)
* **Florian Schmidt, Jakob Zimmermann, Tanmay Tanna et al.** *Noninvasive assessment of gut function using transcriptional recording sentinel cells.* Science, 376:eabm6038, 2022. (https://www.science.org/doi/10.1126/science.abm6038)
* **Tanmay Tanna, Florian Schmidt, Maria Y. Cherepkova et al.** *Recording transcriptional histories using Record-seq.* Nature Protocols, 15:513–539, 2020. (https://www.nature.com/articles/s41596-019-0253-4)

## Repository Organization

The Record-seq software ecosystem is partitioned into three modules:

* **[Primary Analyis Workflow](https://github.com/plattlab/Primary-Analysis-Workflow):** A Snakemake workflow for sample processing and generating count tables.
* **[Secondary Analysis Package](https://github.com/plattlab/recoRdseq):** An R package supporting data processing, differential expression (DE) testing, and visualization of results.
* **[Spacer Acquisition Modelling](https://github.com/plattlab/Spacer-Acquisition-Modelling):** Characterization and modelling of position-specific Record-seq spacer acquisition rates.

## Reproduction of Figures

The `scripts/` directory contains the code required to reproduce the specific computational benchmarks described in the manuscript.

* **Figure 2 (Runtime Performance):** Benchmarks the WFA2-based `PatternBoundExtractor` against previous fuzzysearch methods to demonstrate improved scalability. **Implementation:** `placeholder.py`.
* **Figure 3 (TU-based Feature Counting):** Compares signal recovery and statistical power between transcription unit (TU) and gene-body-centric quantification. **Implementation:** `placeholder.Rmd`.
* **Figure 4 (Acquisition Bias Modelling):** Scripts for normalized profile analysis and model evaluation are hosted within the **[Spacer Acquisition Modelling](https://github.com/plattlab/Spacer-Acquisition-Modelling)** repository.

## Installation for Figure Reconstruction

Use these instructions to establish the environment required to run the figure reconstruction and benchmarking scripts located in the `scripts/` directory.

```bash
# Clone the repository
git clone [https://github.com/plattlab/An-End-To-End-Computational-Framework-for-Record-seq-Transcriptional-Recording-Data.git]
cd An-End-To-End-Computational-Framework-for-Record-seq-Transcriptional-Recording-Data

# Establish the environment for figure reproduction
mamba env create -f envs/env_fig_repro.yaml
conda activate recordseq-fig-repro
```

## Citation

Please cite the following publication when using this framework:

> **"An End-To-End Computational Framework for 'Record-seq' Transcriptional Recording Data"**
> *Journal placeholder* 
