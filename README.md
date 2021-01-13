# Introduction
The main purpose of this project is to investigate the feasibility and effectiveness of comparing cell types across different species using gene expression atlases constructed via PCA. To accomplish the goal, several experiments will be conducted, including: 

- Projecting the cell samples from a species to the gene expression atlas of another species using common orthologs and evaluate the projection by comparing the distances of the correlated cell types between the two species. 
- Analyse the technical and genetic factors that contribute to the similarity and discrepancy of the projection regions of the correlated cell types from the two species.
- Explore the limit of PCA that is able to grasp the major structure of sample distribution by controling the number and the variability of genes included. 

Due to high genome similarity and abundance of genetic data for both human and mouse blood cells, the experiment will compare human and mouse blood cells by projecting mouse samples to the human blood atlas.

# Plan
The data analysis of this project can be roughly divided into three steps: 
1. Construct human BloodAtlas and project test samples from human and mouse onto the BloodAtlas.
2. Exploratory analyses to find factors that might explain the area of projection of the mouse test samples in relation to the human BloodAtlas.
3. Case-control experiments on the gene expression data to find marginal values of key parameters that maintain the major structure of a PCA atlas. 

# Notebooks
The data analysis process conducted in this project will be recorded in a serie of jupyter notebooks with each of them handle a relative smaller analysis task. The index of these notebooks are listed below.

1.1-build-human-BloodAtlas-with-S4M-code
1.2-project-human-test-samples-to-hBloodAtlas
1.3-make-ortholog-map-between-human-and-mouse-ensembleID
1.4-project-mouse-test-samples-to-hBloodAtlas
2.1-compare-distances-between-projected-mcells-and-cells-in-hBloodAtlas
2.2-construct-mouse-BloodAtlas-via-Goodell-data
2.3-compare-expression-variation-between-mouse-and-human-bloodcells

# Datasets
| Name | Source | Note |
| ------ | ------ | ------ |
| Stemformatics Blood Atlas | https://www.stemformatics.org/atlas/blood | Contain samples matrix, genes matrix, expressions matrix & colors. |
| Haemopedia-Human-RNASeq_raw | https://haemosphere.org/datasets/show | Human haematopoietic cell RNASeq expression data from the Hilton Laboratory at the Walter and Eliza Hall Institute. |
| Haemopedia-Mouse-RNASeq_raw | https://haemosphere.org/datasets/show | Mouse RNASeq Atlas created by the Hilton lab at the Walter and Eliza Hall Institute on sorted wildtype mouse haematopoietic cells. |
| Goodell Mouse dataset | https://haemosphere.org/datasets/show | This is a collection of microarray expression data that covers 10 different murine haematopoietic cell populations which have been FACS sorted. It has been generated by Peggy Goodell's lab in Baylor College of Medicine. |
| HOM_MouseHumanSequence | http://www.informatics.jax.org/downloads/reports/ | Human and Mouse Homology Classes with Sequence information (tab-delimited) |