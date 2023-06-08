# Master 2 in Bioinformatics, University of Rennes, Internship project.

Important note : **Not all scricpt has been cleaned and properly commented yet**

Script used in my 2023 Master degree internship project.

This aimed to build co-expression network from single cell RNA-seq data.

Input to the co-expression computation is a R `list` with 2 elements :
 * (named ctmat) A read counts matrix with cells in column (named) and genes in rows (named).
 * (neamed meta) A data table with a row for each cells and column with : cells labels (corresponding to count matrix columns names), cell-types, subject labels. columns names have to respectevely be cell, cell_type, and subject.

scripts used to made the figures of the internship report related to that project can be found in `figures`

Data used can be found at the Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/) (GEO) archives, Synapse archive(https://www.synapse.org/#), or through the links provided in the table below. Noticing that ROSMAP data had to be requested and that metadata files related to the Lau dataset were requested to the authors of the related publication and were not available in the GEO archive at the date of the study. The GO annotations file used for the GBA analysis is the human annotation file of the 2022-12-04 release and can be found in the GO FTP archive (http://release.geneontology.org/).

| Dataset name | GEO/SynID/link                    |
| ------------ | --------------------------------- |
| Lau          | GSE157827                         |
| Lim          | GSE180928                         |
| Nagy         | GSE144136                         |
| Pineda       | GSE174332                         |
| Ramos        | GSE217511                         |
| ROSMAP       | syn18485175                       |
| Velmeshev    | https://cells.ucsc.edu/?ds=autism |
