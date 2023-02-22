FROM rocker/shiny-verse:latest
MAINTAINER Kan Etoh
RUN apt-get update && apt-get install -y \
    build-essential \
    libglpk40 \
    libbz2-dev \
    liblzma-dev \
    libgsl-dev \
    r-cran-gsl \
    wget \
    git

RUN R -e "install.packages('BiocManager',repos='http://cran.rstudio.com/')" && \
    R -e "BiocManager::install('devtools', update = F)" && \
    R -e "devtools::install_github('haizi-zh/bedtorch')" && \
    R -e "devtools::install_github('YuLab-SMU/clusterProfiler.dplyr')" && \
    R -e "BiocManager::install('shiny', update = F)" && \
    R -e "BiocManager::install('DT', update = F)" && \
    R -e "BiocManager::install('gdata', update = F)" && \
    R -e "BiocManager::install('rstatix', update = F)" && \
    R -e "BiocManager::install('multcomp', update = F)" && \
    R -e "BiocManager::install('ggrepel', update = F)" && \
    R -e "BiocManager::install('ggdendro', update = F)" && \
    R -e "BiocManager::install('ggplotify', update = F)" && \
    R -e "BiocManager::install('gridExtra', update = F)" && \
    R -e "BiocManager::install('cowplot', update = F)" && \
    R -e "BiocManager::install('DESeq2', update = F)" && \
    R -e "BiocManager::install('ggnewscale', update = F)" && \
    R -e "BiocManager::install('org.Hs.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Mm.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Rn.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Dm.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Ce.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Xl.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Bt.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Cf.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Dr.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Gg.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Mmu.eg.db', update = F)" && \
    R -e "BiocManager::install('org.Pt.eg.db', update = F)" && \
    R -e "BiocManager::install('biomaRt', update = F)" && \
    R -e "BiocManager::install('DOSE', update = F)" && \
    R -e "BiocManager::install('msigdbr', update = F)" && \
    R -e "BiocManager::install('genefilter', update = F)" && \
    R -e "BiocManager::install('ComplexHeatmap', update = F)" && \
    R -e "BiocManager::install('shinyBS', update = F)" && \
    R -e "BiocManager::install('plotly', update = F)" && \
    R -e "BiocManager::install('shinyjs', update = F)" && \
    R -e "BiocManager::install('dorothea', update = F)" && \
    R -e "BiocManager::install('ggpubr', update = F)" && \
    R -e "BiocManager::install('enrichplot', update = F)" && \
    R -e "BiocManager::install('clusterProfiler', update = F)" && \
    R -e "BiocManager::install('GenomicRanges', update = F)" && \
    R -e "BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Dmelanogaster.UCSC.dm6.ensGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Rnorvegicus.UCSC.rn6.refGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Celegans.UCSC.ce11.refGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Btaurus.UCSC.bosTau8.refGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Cfamiliaris.UCSC.canFam3.refGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Drerio.UCSC.danRer10.refGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Ggallus.UCSC.galGal4.refGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Mmulatta.UCSC.rheMac8.refGene', update = F)" && \
    R -e "BiocManager::install('TxDb.Ptroglodytes.UCSC.panTro4.refGene', update = F)" && \
    R -e "BiocManager::install('rtracklayer', update = F)" && \
    R -e "BiocManager::install('ChIPseeker', update = F)" && \
    R -e "BiocManager::install('ChIPpeakAnno', update = F)" && \
    R -e "BiocManager::install('rGREAT', update = F)" && \
    R -e "BiocManager::install('tidyverse', update = F)" && \
    R -e "BiocManager::install('Gviz', update = F)" && \
    R -e "BiocManager::install('FindIT2', update = F)" && \
    R -e "BiocManager::install('limma', update = F)" && \
    R -e "BiocManager::install('ggseqlogo', update = F)" && \
    R -e "BiocManager::install('marge', update = F)" && \
    R -e "BiocManager::install('Rsubread', update = F)" && \
    R -e "devtools::install_github('ColeWunderlich/soGGi')" && \
    R -e "devtools::install_github('robertamezquita/marge',ref = 'master')" && \
    R -e "devtools::install_github('amitjavilaventura/plotmics')" && \
    R -e "BiocManager::install('colorspace', update = F)" && \
    R -e "BiocManager::install('ggcorrplot', update = F)" && \
    R -e "BiocManager::install('RColorBrewer', update = F)" && \
    R -e "BiocManager::install('plyranges', update = F)" && \
    R -e "BiocManager::install('ggvenn', update = F)"
##Remove the unnecessary genomes for HOMER
RUN mkdir -p /srv/shiny-server/EpigenomeChef && \
    rm -rf /srv/shiny-server/hello && \
    mkdir -p /usr/local/homer && \
    cd /usr/local/homer && \
    wget http://homer.ucsd.edu/homer/configureHomer.pl && \
    perl configureHomer.pl -install && \
    perl configureHomer.pl -install hg19 && \
    perl configureHomer.pl -install hg38 && \
    perl configureHomer.pl -install mm10 && \
    perl configureHomer.pl -install rn6 && \
    perl configureHomer.pl -install ce11 && \
    perl configureHomer.pl -install dm6 && \
    perl configureHomer.pl -install canFam3 && \
    perl configureHomer.pl -install galGal4 && \
    perl configureHomer.pl -install danRer10 && \
    perl configureHomer.pl -install rheMac8 && \
    perl configureHomer.pl -install panTro4
COPY ui.R /srv/shiny-server/EpigenomeChef/
COPY server.R /srv/shiny-server/EpigenomeChef/
COPY global.R /srv/shiny-server/EpigenomeChef/
COPY navAppend.js /srv/shiny-server/EpigenomeChef/
COPY navAppend2.js /srv/shiny-server/EpigenomeChef/
COPY data /srv/shiny-server/EpigenomeChef/data/
COPY www /srv/shiny-server/EpigenomeChef/www/
COPY shiny-server.conf /etc/shiny-server/
ENV PATH $PATH:/usr/local/homer/bin
ENV PERL5LIB $PERL5LIB:/usr/local/homer/bin
RUN chown -R shiny:shiny /srv/shiny-server && \
    chown -R shiny:shiny /usr/local/homer
EXPOSE 3838
CMD exec shiny-server >> /var/log/shiny-server.log 2>&1
