# Transformer-based Analysis of Diabetic Foot Ulcer (DFU) Single-cell RNA-seq Data

This project utilizes a transformer model to analyze a single-cell RNA-seq dataset from human foot skin, distinguishing between Healthy, Diabetic, Healing DFU, and Non-Healing DFU conditions. The goal is to identify the most important gene features for each condition using attention scores generated from the transformer model.

## Methods

### Data
The analysis is based on the single-cell RNA-seq dataset from human foot skin samples, as published by [Theocharidis et al., 2022](https://doi.org/10.1038/s41467-021-27884-3). We analyzed:
- **Control samples**: Patients without comorbidities (4 samples)
- **Diabetic samples**: 8 diabetic samples without ulcer presence
- **Healing DFU**: 8 samples with healing DFU
- **Non-Healing DFU**: 4 samples with non-healing DFU

### Preprocessing
Each sample underwent individual preprocessing:
1. **Cell Filtering**: Cells with fewer than 200 genes and genes expressed in fewer than 3 cells were removed.
2. **Count and Mitochondrial Filtering**: Only cells with n_Counts < 25,000 and mitochondrial percentage < 15% were retained.
3. **Integration**: Samples were integrated using Seurat V4 [Hao et al., 2021]. 
4. **Normalization and Feature Selection**: SCTransform was applied for normalization, scaling, and selection of highly variable genes (989 genes).
5. **Dimensionality Reduction and Clustering**: PCA was performed, retaining 60 principal components for UMAP visualization and clustering via the Louvain algorithm.

### Cell Typing
After clustering, cell types were identified based on known markers and annotations.

## Transformer Model
To classify between conditions, a transformer model was implemented following the structure outlined by [Zhang et al., 2022].

- **Architecture**: Three attention layers with three attention heads each, taking in the normalized gene expression of each cell.
- **Self-Attention Mechanism**: For each gene `g` in a cell `i`, the query, key, and value vectors were computed using learned weights, representing the gene interactions within the cell.
- **Attention Scores Calculation**: Attention weights were calculated between genes using a softmax layer, obtaining a probability distribution that reflects gene-gene interactions.
- **Representation of Query Genes**: The attention outputs for each gene were combined into a final set of representations, with concatenated embeddings for additional factors such as Age, Gender, BMI, Creatinine, and Urea.
- **Training**: The model was trained for 200 epochs with early stopping, a dropout rate of 0.4, and an Adam optimizer with learning rates between 0.001 and 0.0001 depending on cell type.

### Model Training
- **Dataset Split**: 70% of the data was used for training, 10% for validation, and the remaining for testing.
- **Optimizer**: Adam with weight decay (0.01) and negative log-likelihood loss function.
- **Batch Size**: 16.

### Pathway Enrichment Analysis
Attention scores from the final layer were summed for each gene and sorted to identify the top 20 genes. Enrichment analysis was conducted using the **enrichr** package, with pathway databases including Hallmark, Gene Ontology, KEGG, Reactome, and Diabetes Perturbations GEO.

## Results
Attention weights allowed us to identify significant gene pathways and features for each condition, contributing to a better understanding of molecular functions related to DFU healing.

## References
1. Theocharidis et al., 2022 - [Single cell transcriptomic landscape of diabetic foot ulcers](https://doi.org/10.1038/s41467-021-27884-3)
2. Hao et al., 2021 - "Integrated analysis of multimodal single-cell data"
3. Hafemeister and Satija, 2019 - "Normalization and variance stabilization of single-cell RNA-seq data"
4. Zhang et al., 2022 - "Transformer for gene expression modeling (T-GEM)"
5. Ba et al., 2016 - "Layer normalization"
6. Kuleshov et al., 2016 - "Enrichr: a comprehensive gene set enrichment analysis web server"
7. Liberzon et al., 2015 - "The molecular signatures database hallmark gene set collection"
8. Aleksander et al., 2023 - "The gene ontology knowledgebase"
9. Kanehisa and Goto, 2000 - "KEGG: kyoto encyclopedia of genes and genomes"
10. Croft et al., 2010 - "Reactome: a database of reactions, pathways and biological processes"

---

This repository provides a foundation for applying transformers in the analysis of single-cell RNA-seq data, highlighting feature selection and pathway enrichment for diabetic foot ulcer conditions.

