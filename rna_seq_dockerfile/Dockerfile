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

# Install bioconductor packages for microarray analysis
RUN Rscript -e "BiocManager::install('oligo')"
RUN Rscript -e "BiocManager::install('hta20transcriptcluster.db')"
RUN Rscript -e "BiocManager::install('affycoretools')"
RUN Rscript -e "BiocManager::install('pd.hta.2.0')"

#Install HTSeq-count
RUN pip install HTSeq

RUN Rscript -e "install.packages('VennDiagram')"
RUN Rscript -e "install.packages('ggrepel')"

#Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar xvjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && ./configure --without-curses --disable-bz2 --disable-lzma \
    && make \
    && make install


# Install Ballgown
RUN Rscript -e "BiocManager::install('ballgown')"

