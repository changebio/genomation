# RRBS from ENCODE
rrbs.IMR90.path <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/wgEncodeHaibMethylRrbsImr90UwSitesRep1.bed.gz"
rrbs.H1.path <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/wgEncodeHaibMethylRrbsH1hescHaibSitesRep1.bed.gz"

# read the RRBS data to GRanges
library(genomation)

rrbs.IMR90<-readGeneric(rrbs.IMR90.path, chr = 1, start = 2, end = 3, strand = 6,meta.cols = list(score=11), keep.all.metadata = FALSE, zero.based = TRUE, skip = 1)
rrbs.H1<-readGeneric(rrbs.H1.path, chr = 1, start = 2, end = 3, strand = 6,meta.cols = list(score=11), keep.all.metadata = FALSE, zero.based = TRUE, skip = 1)

# mapping to [0, 1] range
rrbs.IMR90$score<-rrbs.IMR90$score/100
rrbs.H1$score<-rrbs.H1$score/100


# download WGBS data from Roadmap Epigenomics, convert to GRanges
library(rtracklayer)
wgbs.H1.path <- "http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/E003_WGBS_FractionalMethylation.bigwig"
wgbs.IMR90.path <- "http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/E017_WGBS_FractionalMethylation.bigwig"
grange <- GRanges("chr3", IRanges(1, 2e8)) #chr3
wgbs.H1 <- import.bw(wgbs.H1.path, selection=BigWigSelection(grange))
wgbs.IMR90 <- import.bw(wgbs.IMR90.path, selection=BigWigSelection(grange))



# Then we read transcripts from the refseq bed file 
# within 4 bp distance from the longest transcript of each gene.
refseq.path <- "/mnt/local-disk1/rsgeno2/MAmotif/RACK7/routput/chr3.refseq.hg19.bed"
transcriptFeat=readTranscriptFeatures(refseq.path, up.flank = 2000, down.flank = 2000)

# Once we have read the files
# now we can build base-pair resolution matrices of scores for each experiment.
# The returned list of matrices can be used to draw meta profiles or heatmaps of read coverage around promoters.


targets <- list(WGBS.H1=wgbs.H1, WGBS.IMR90=wgbs.IMR90, RRBS.H1=rrbs.H1, RRBS.IMR90=rrbs.IMR90) 
sml <- ScoreMatrixList(targets=targets, windows=transcriptFeat$promoters, bin.num=20, strand.aware=TRUE, weight.col="score", is.noCovNA = TRUE)

sml.sub =intersectScoreMatrixList(sml, reorder = FALSE) # because scoreMatrices have different dimensions, we need them to cover the same promoters for the heatmap
multiHeatMatrix(sml.sub, xcoords=c(-2000, 2000), matrix.main=names(targets), xlab="region around TSS")
