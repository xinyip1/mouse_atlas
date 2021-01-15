# Introduction
The main purpose of this project is to build a prototype mouse blood atlas to serve as a foundation of estabilishing a comprehensive reference map for the mouse blood development. As well as, aid with data reusability by making the atlas compatible with data derived from various technical platforms. The goal is approached by following steps:

1. Collect relevant transcriptional data from multiple accessible public databases, including Stemformatics(https://www.stemformatics.org), Encode(https://www.encodeproject.org), and Haemosphere (https://haemosphere.org) data portals.
2. Integrat the collected gene expression data and metadata, including sample reannotation to increase consistency of metadata.
3. Conduct variance analysis to access variance dependence of each gene due to platform, and control for batch effect through gene filtering.
4. Visulize the processed data via dimension reduction method PCA.

Detailed code implementation the sequential data processing is stored in the /notebook folder.

# Data
Data used in this project is stored under /data folder, in which they are further divided into three categories: 

- raw: datasets that are downloaded from external websites
- interim: intermediate data that has been transformed
- processed: the final datasets for modelling 

