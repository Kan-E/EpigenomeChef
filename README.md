# EpigenomeChef
EpigenomeChef is a platform of systematic epigenome data analysis which can automatically detect, integrate, and visualize the epigenetic information and its interaction with transcriptome without bioinformatics skills.<br>
Data downloaded from ChIP-atlas, a database of ChIP-seq, ATAC-seq, and Bisulfite-seq, can be used as input.<br>
In addition, integrated analysis with transcriptome data is possible by uploading the results file of differentially expressed gene (DEG) analysis obtained from RNA-seq data.<br>

# Local installation
EpigenomeChef can be used by installing the Docker image.<br>

0. (If you do not have a Docker environment) Install Docker <br>
1. Run the following command (Once you run it, you won't need it again) <br> 
```
#AMD64 architecture (Windows, Linux, MacOS(Intel))
##HOMER human (hg19) 
docker pull kanetoh1030/shiny-epigenomechef:hg19v0.0.1-amd64
##HOMER mouse (mm10) 
docker pull kanetoh1030/shiny-epigenomechef:mm10v0.0.1-amd64

#ARM64 architecture (MacOS (M1/M2))
##HOMER human (hg19) 
docker pull kanetoh1030/shiny-epigenomechef:hg19v0.0.1-arm64
##HOMER mouse (mm10) 
docker pull kanetoh1030/shiny-epigenomechef:mm10v0.0.1-arm64

```
2. Run the following command for the launch EpigenomeChef on the browser
```
#AMD64 architecture (Windows, Linux, MacOS(Intel))
##HOMER human (hg19)
docker run --rm -p 3838:3838 kanetoh1030/shiny-epigenomechef:hg19v0.0.1-amd64
##HOMER mouse (mm10) 
docker run --rm -p 3838:3838 kanetoh1030/shiny-epigenomechef:mm10v0.0.1-amd64

#ARM64 architecture (MacOS (M1/M2))
##HOMER human (hg19) 
docker run --rm -p 3838:3838 kanetoh1030/shiny-epigenomechef:hg19v0.0.1-arm64
##HOMER mouse (mm10) 
docker run --rm -p 3838:3838 kanetoh1030/shiny-epigenomechef:mm10v0.0.1-arm64

```
Access [http://localhost:3838](http://localhost:3838).

If you need help, please create an issue on [Github](https://github.com/Kan-E/EpigenomeChef/issues) or [contact us](mailto:kaneto@kumamoto-u.ac.jp). <br>

# Reference
Shiny framework
- Winston Chang, Joe Cheng, JJ Allaire, Carson Sievert, Barret Schloerke, Yihui
  Xie, Jeff Allen, Jonathan McPherson, Alan Dipert and Barbara Borges (2021).
  shiny: Web Application Framework for R. R package version 1.7.1.
  https://CRAN.R-project.org/package=shiny
- Eric Bailey (2022). shinyBS: Twitter Bootstrap Components for Shiny. R package
  version 0.61.1. https://CRAN.R-project.org/package=shinyBS
- Yihui Xie, Joe Cheng and Xianying Tan (2022). DT: A Wrapper of the JavaScript
  Library 'DataTables'. R package version 0.23.
  https://CRAN.R-project.org/package=DT
  
DESeq2 and limma (for differential accessibility region analysis)
- Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for
  RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
- Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses for RNA-sequencing
  and microarray studies. Nucleic Acids Research 43(7), e47.

ggdendro (for dendrograms)
- Andrie de Vries and Brian D. Ripley (2020). ggdendro: Create Dendrograms and Tree Diagrams Using 'ggplot2'. R package version 0.1.22. https://CRAN.R-project.org/package=ggdendro

ggcorrplot (for correlation plot)
- Kassambara A (2022). _ggcorrplot: Visualization of a Correlation Matrix using 'ggplot2'_. R package version 0.1.4,
  <https://CRAN.R-project.org/package=ggcorrplot>.

rGREAT, clusterProfiler, DOSE, msigdbr, dorothea (for enrichment analysis)
- Gu Z (2022). _rGREAT: GREAT Analysis - Functional Enrichment on Genomic Regions_. https://github.com/jokergoo/rGREAT,
  http://great.stanford.edu/public/html/.
- T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
- Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015 31(4):608-609
- Dolgalev I (2022). _msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format_. R
  package version 7.5.1, <https://CRAN.R-project.org/package=msigdbr>.
- Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. 'Benchmark and integration of resources for the estimation of human transcription factor activities.' Genome Research. 2019. DOI: 10.1101/gr.240663.118.

AnnotationDbi, org.Hs.eg.db, org.Mm.eg.db (for genome wide annotation)
- Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2020). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.52.0. https://bioconductor.org/packages/AnnotationDbi
- Marc Carlson (2020). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.12.0.
- Marc Carlson (2020). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.12.0.

genefilter (for z-score normalization)
- R. Gentleman, V. Carey, W. Huber and F. Hahne (2021). genefilter: methods for filtering genes from high-throughput experiments. R package version 1.72.1.

ComplexHeatmap (for heatmap and k-means clustering)
- Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.

ggplot2 and ggpubr (for boxplot and scater plot)
- H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
- Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.4.0. https://CRAN.R-project.org/package=ggpubr

venn (for venn diagram analysis)
- Mitjavila-Ventura A (2022). _plotmics: Visualization of Omics And Sequencing Data In R_. https://amitjavilaventura.github.io/plotmics/index.html,
  https://github.com/amitjavilaventura/plotmics.

GenomicRanges, soGGi, ChIPQC,ChIPseeker,TxDb.Mmusculus.UCSC.mm10.knownGene, TxDb.Hsapiens.UCSC.hg19.knownGene (for manupilate genomic information and annotation)
- Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118. doi:10.1371/journal.pcbi.1003118"
- Dharmalingam G, Carroll T (2015). _soGGi: Visualise ChIP-seq, MNase-seq and motif occurrence as aggregate plots Summarised Over Grouped Genomic Intervals_. R package version 1.10.4.
- Thomas Samuel Carroll, Ziwei Liang, Rafik Salama, Rory Stark, Ines de Santiago: Impact of artefact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data. Frontiers in Genetics, in press.
- Guangchuang Yu, Li-Gen Wang, and Qing-Yu He. ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. Bioinformatics 2015, 31(14):2382-2383
- Lihua J Zhu, Claude Gazin, Nathan D Lawson, Herve Pages, Simon M Lin, David S Lapointe and Michael R Green, ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics. 2010, 11:237
- Team BC, Maintainer BP (2019). _TxDb.Mmusculus.UCSC.mm10.knownGene: Annotation package for TxDb object(s)_. R package version 3.10.0.
- Carlson M, Maintainer BP (2015). _TxDb.Hsapiens.UCSC.hg19.knownGene: Annotation package for TxDb object(s)_. R package version 3.2.2.
- Team BC, Maintainer BP (2022). _TxDb.Hsapiens.UCSC.hg38.knownGene: Annotation package for TxDb object(s)_. R package version 3.15.0.

dplyr and tidyr (for data manipulation)
- Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr
- Hadley Wickham (2021). tidyr: Tidy Messy Data. R package version 1.1.3. https://CRAN.R-project.org/package=tidyr

rtracklayer, Rsubread, Rsamtools
- M. Lawrence, R. Gentleman, V. Carey: "rtracklayer: an {R} package for interfacing with genome browsers". Bioinformatics 25:1841-1842.
- Liao Y, Smyth GK and Shi W (2019). The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads.
  Nucleic Acids Research 47(8), e47.
- Morgan M, Pagès H, Obenchain V, Hayden N (2022). _Rsamtools: Binary alignment (BAM), FASTA, variant call (BCF), and tabix file import_. R package version
  2.12.0, <https://bioconductor.org/packages/Rsamtools>.

FindIT2 (for with RNAseq)
- Shang, G.-D., Xu, Z.-G., Wan, M.-C., Wang, F.-X. & Wang, J.-W.  FindIT2: an R/Bioconductor package to identify influential transcription factor and
  targets based on multi-omics data.  BMC Genomics 23, 272 (2022)

marge, ggseqlogo (for motif analysis)
- Amezquita R (2022). _marge: API for HOMER in R for Genomic Analysis using Tidy Conventions_. R package version 0.0.4.9999.
- Wagih O (2017). _ggseqlogo: A 'ggplot2' Extension for Drawing Publication-Ready Sequence Logos_. R package version 0.1,
  <https://CRAN.R-project.org/package=ggseqlogo>.
  
colorspace
- Zeileis A, Hornik K, Murrell P (2009). "Escaping RGBland: Selecting Colors for Statistical Graphics." _Computational Statistics \& Data Analysis_,
  *53*(9), 3259-3270. doi:10.1016/j.csda.2008.11.033 <https://doi.org/10.1016/j.csda.2008.11.033>.
  
bedtorch (for bedtools)
- Zheng H (2021). _bedtorch: R package for fast BED-file manipulation_. R package version 0.1.12.12, <https://github.com/haizi-zh/bedtorch>.

Gviz (for track plot)
- Hahne F, Ivanek R. Visualizing Genomic Data Using Gviz and Bioconductor. Methods Mol Biol. 1418:335-51 (2016).

# License
This shiny code is licensed under the GPLv3. Please see the file [LICENSE.md](https://github.com/Kan-E/EpigenomeChef/blob/main/LICENSE.md) for information.<br>
```
EpigenomeChef
Shiny App for automated, systematic, and integrated epigenetic analysis
Copyright (C) 2022  Kan Etoh

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

You may contact the author of this code, Kan Etoh, at <kaneto@kumamoto-u.ac.jp>
```
# Author

Kan Etoh
<[kaneto@kumamoto-u.ac.jp](mailto:kaneto@kumamoto-u.ac.jp)>
