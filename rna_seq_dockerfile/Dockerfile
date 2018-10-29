#Use Bioconductor release_core2 as the base
#Version R3.5.1_Bioc3.7 is the latest as of 29/10/2018
from bioconductor/release_core2:R3.5.1_Bioc3.7

#Install Bioconductor packages necessary for analysis
RUN Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('DESeq2')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('tximport')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R'); biocLite('EnsDb.Hsapiens.v86')"

#Install non-bioconductor packages necessary for analysis
RUN Rscript -e "install.packages('tidyverse')"
RUN Rscript -e "install.packages('RColorBrewer')"
RUN Rscript -e "install.packages('pheatmap')"

#Install HTSeq-count
RUN pip install HTSeq

