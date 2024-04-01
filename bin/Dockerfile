FROM bioconductor/bioconductor_docker:RELEASE_3_17

LABEL maintainer="Victoria Flanary <flanary@uab.edu>"
LABEL description="An environment for analyzing ATAC-seq data" 

### CRAN Packages
RUN R -e 'install.packages(c("tidyverse", "plyr", "remotes",  "BiocManager", "devtools", "lintr", "styler", "here", "markdown", "rmarkdown"))'
RUN R -e 'install.packages(c("ggthemes", "ggpubr", "patchwork", "colorspace", "RColorBrewer"))'
RUN R -e 'install.packages(c("gprofiler2", "Matrix"))'
### BioConductor Packages
#### Needed to complete exercises for Computational Genomics in R
RUN R -e 'BiocManager::install(c("qvalue", "plot3D", "pheatmap", "ComplexHeatmap", "cowplot"))'
RUN R -e 'BiocManager::install(c("cluster", "NBClust", "fastICA", "NMF", "matrixStats"))'
RUN R -e 'BiocManager::install(c("Rtsne","mosaic", "knitr", "genomation", "BiocGenerics"))'
RUN R -e 'BiocManager::install(c("ggbio","Gviz", "DESeq2", "RUVSeq", "gProfileR", "clusterProfiler"))'
RUN R -e 'BiocManager::install(c("ggfortify", "corrplot", "gage", "EDASeq", "citr", "formatR"))'
RUN R -e 'BiocManager::install(c("svglite", "Rqc", "ShortRead", "QuasR", "methylKit"))'
RUN R -e 'BiocManager::install(c("FactoMineR", "iClusterPlus", "enrichR", "caret", "xgboost"))'
RUN R -e 'BiocManager::install(c("glmnet", "DALEX", "kernlab", "pROC", "nnet", "RANN", "ranger"))'
RUN R -e 'BiocManager::install(c("GenomeInfoDb", "GenomicRanges", "GenomicAlignments", "ChIPQC"))'
RUN R -e 'BiocManager::install(c("circlize", "rtracklayer", "ChIPseeker", "DiffBind"))'
RUN R -e 'BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg19"))'
RUN R -e 'BiocManager::install(c("tidyr", "AnnotationHub", "GenomicFeatures", "normr"))'
RUN R -e 'BiocManager::install(c("MotifDb", "TFBSTools", "rGADEM", "JASPAR2018"))'
RUN R -e 'BiocManager::install(c("ensembldb", "AnnotationDbi", "org.Hs.eg.db"))'

### GitHub Packages
RUN R -e 'devtools::install_github("demuellae/muLogR")'
RUN R -e 'devtools::install_github("demuellae/muRtools")'
RUN R -e 'devtools::install_github("demuellae/muReportR")'
RUN R -e 'devtools::install_github("GreenleafLab/ChrAccR", dependencies=TRUE)'