# Introduction
The main purpose of this project is to build a prototype mouse blood atlas to serve as a foundation of estabilishing a comprehensive reference map for the mouse blood development. As well as, aid with data reusability by making the atlas compatible with data derived from various technical platforms. The goal is approached by following steps:

1. Collect relevant transcriptional data from multiple accessible public databases, including Stemformatics(https://www.stemformatics.org), Encode(https://www.encodeproject.org), and Haemosphere (https://haemosphere.org) data portals.
2. Integrat the collected gene expression data and metadata, including sample reannotation to increase consistency of metadata.
3. Conduct variance analysis to access variance dependence of each gene due to platform, and control for batch effect through gene filtering.
4. Visulize the processed data via dimension reduction method PCA.

# Notebooks
The data analysis process conducted in this project are recorded in a serie of jupyter notebooks in /notebook folder, with each of them handle a relative smaller analysis task. The index of these notebooks are listed below.

- 1_data_collection_Encode
- 1_data_collection_Haemosphere
- 1_data_collection_Stemformatics
- 2_data_integration
- 3_threshold_selection_and_data_integration
- 4_PCA
- 5_test_stability

# Data Folder Structure
The data used in this project is stored under /data folder, in which they are further divided into three categories:
- raw: datasets that are directly downloaded from external websites
- interim: intermediate data that has been transformed 
- processed: the final datasets for modelling

# Datasets 
Details of the included datasets are listed below. 
| Dataset Name            | GEO/ArrayExpression Accession & Data Link                                  | Platform                         |
| ----------------------- | -------------------------------------------------------------------------- | -------------------------------- |
| Copley                  | [GSE41758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41758)    | Affymetrix Mouse Gene 1 ST Array |
| Immgen                  | [GSE15907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907)    | Affymetrix Mouse Gene 1 ST Array |
| Mackay                  | [GSE47045](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47045)    | Affymetrix Mouse Gene 1 ST Array |
| Pandya                  | [GSE47605](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47605)    | Affymetrix Mouse Gene 1 ST Array |
| Polo                    | [GSE22043](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22043)    | Affymetrix Mouse Gene 1 ST Array |
| Zigmond                 | [GSE55606](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55606)    | Affymetrix Mouse Gene 1 ST Array |
| Anandasabapathy         | [GSE29949](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29949)    | Affymetrix Mouse430 2 Array      |
| Beerman                 | [GSE55525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55525)    | Affymetrix Mouse430 2 Array      |
| Gene-Expression-Commons | [GSE34723](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34723)    | Affymetrix Mouse430 2 Array      |
| Goodell                 | [GSE6506](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6506)      | Affymetrix Mouse430 2 Array      |
| Martinez                | [GSE35435](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35435)    | Affymetrix Mouse430 2 Array      |
| Mikkelsen               | [GSE10871](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10871)    | Affymetrix Mouse430 2 Array      |
| Rossi                   | [GSE4332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4332)      | Affymetrix Mouse430 2 Array      |
| Cabezas-Wallscheid      | [E-MTAB-4549](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-4549/) | Illumina HiSeq                   |
| Haemopedia_RNA-seq      | [GSE116177](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116177)  | Illumina HiSeq                   |
| Takata                  | [GSE99078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99078)    | Illumina HiSeq                   |
| Chaudhury               | [GSE107662](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107662)  | Illumina NextSeq                 |
| Immgen_ULI              | [GSE109125](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109125)  | Illumina NextSeq                 |
| ENCODE                  | [Link](https://www.encodeproject.org/search/?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&files.file_type=tsv&audit.WARNING.category!=low+replicate+concordance&status=released&assay_slims=Transcription&biosample_ontology.classification!=cell+line&biosample_ontology.cell_slims=leukocyte&biosample_ontology.cell_slims=hematopoietic+cell&biosample_ontology.cell_slims=T+cell&biosample_ontology.cell_slims=myeloid+cell&biosample_ontology.cell_slims=monocyte&biosample_ontology.cell_slims=stem+cell&biosample_ontology.cell_slims=progenitor+cell&award.project=ENCODE)                                  | Illumina NextSeq/HiSeq            |              
| GGR                     | [Link](https://www.encodeproject.org/search/?type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&files.file_type=tsv&audit.WARNING.category%21=low+replicate+concordance&status=released&assay_slims=Transcription&biosample_ontology.classification%21=cell+line&biosample_ontology.cell_slims=leukocyte&biosample_ontology.cell_slims=hematopoietic+cell&biosample_ontology.cell_slims=T+cell&biosample_ontology.cell_slims=myeloid+cell&biosample_ontology.cell_slims=monocyte&biosample_ontology.cell_slims=stem+cell&biosample_ontology.cell_slims=progenitor+cell&award.project=GGR)                                 | Illumina NextSeq/HiSeq            |        
| Beutner                 | [GSE43808](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43808)    | Illumina MouseWG-6 v2.0           |
| Haemopedia_Plus         | [GSE77098](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77098)    | Illumina MouseWG-6 v2.0           |
| Haniffa                 | [GSE35458](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35458)    | Illumina MouseWG-6 v2.0           |
| Pelekanos               | [GSE31738](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31738)    | Illumina MouseWG-6 v2.0           |


