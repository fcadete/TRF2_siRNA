#Use Bioconductor release_core2 as the base
#Version R3.5.1_Bioc3.7 is the latest as of 29/10/2018
from bioconductor/release_core2:R3.5.1_Bioc3.7

#Install Bioconductor packages necessary for analysis
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('DESeq2')"
RUN Rscript -e "BiocManager::install('tximport')"
RUN Rscript -e "BiocManager::install('EnsDb.Hsapiens.v86')"
RUN Rscript -e "BiocManager::install('org.Hs.eg.db')"
RUN Rscript -e "BiocManager::install('Repitools')"
RUN apt-get install -y libpng-dev
RUN Rscript -e "BiocManager::install('KEGGprofile')"
RUN Rscript -e "BiocManager::install('topGO')"

#Install non-bioconductor packages necessary for analysis
RUN Rscript -e "install.packages('tidyverse')"
RUN Rscript -e "install.packages('RColorBrewer')"
RUN Rscript -e "install.packages('pheatmap')"


#Install HTSeq-count
RUN pip install HTSeq

