library(GenomicRanges)
library(rtracklayer)

reads <- import.bed(con="/mnt/local-disk1/rsgeno2/MAmotif/3.Histone_Broad_hg19/H3K4me3/H1hESC/wgEncodeBroadHistoneH1hescH3k4me3StdAlnRep1.bed", asRangedData=F)

cov <- coverage(reads)

# for a peak on chr1 from 2000:3000
plot(cov[["chr1"]][1000:4000])

library(genomation)
library(genomationData)
genomationDataPath = system.file('extdata',package='genomationData')
bam.files = list.files(genomationDataPath, full.names=TRUE, pattern='bam$')
bam.files = bam.files[!grepl('Cage', bam.files)]

library(GenomicRanges)

ctcf.peaks = readBroadPeak(file.path(genomationDataPath, 'wgEncodeBroadHistoneH1hescCtcfStdPk.broadPeak.gz'))
ctcf.peaks = ctcf.peaks[seqnames(ctcf.peaks) == 'chr21']
ctcf.peaks = ctcf.peaks[order(-ctcf.peaks$signalValue)]
ctcf.peaks = resize(ctcf.peaks, width=1000, fix='center')

sml = ScoreMatrixList(bam.files, ctcf.peaks, bin.num=50, type='bam', cores=2)

# descriptions of file that contain info. about transcription factors
sampleInfo = read.table(system.file('extdata/SamplesInfo.txt', package='genomationData'),header=TRUE, sep='\t')
names(sml) = sampleInfo$sampleName[match(names(sml),sampleInfo$fileName)]

sml1 = sml * 100
sml1
sml[[6]] = sml[[1]]
sml 
sml[[6]] <- NULL
sml.scaled = scaleScoreMatrixList(sml)
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500))

# k-means algorithm with 2 clusters
cl1 <- function(x) kmeans(x, centers=2)$cluster
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl1)

# hierarchical clustering with Ward's method for agglomeration into 2 clusters
cl2 <- function(x) cutree(hclust(dist(x), method="ward"), k=2)
multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl2)

multiHeatMatrix(sml.scaled, xcoords=c(-500, 500), clustfun = cl1, clust.matrix = 1)

plotMeta(mat=sml.scaled, profile.names=names(sml.scaled),xcoords=c(-500, 500),winsorize=c(0,99),centralTend="mean")

plotMeta(mat=sml.scaled, profile.names=names(sml.scaled),xcoords=c(-500, 500),winsorize=c(0,99),centralTend="mean",  smoothfun=function(x) stats::smooth.spline(x, spar=0.5))

plotMeta(mat=sml, profile.names=names(sml),xcoords=c(-500, 500),winsorize=c(0,99),centralTend="mean",  smoothfun=function(x) stats::smooth.spline(x, spar=0.5), dispersion="se", lwd=4)

#ctcf motif from the JASPAR database
ctcf.pfm = matrix(as.integer(c(87,167,281,56,8,744,40,107,851,5,333,54,12,56,104,372,82,117,402, 
                               291,145,49,800,903,13,528,433,11,0,3,12,0,8,733,13,482,322,181, 
                               76,414,449,21,0,65,334,48,32,903,566,504,890,775,5,507,307,73,266, 
                               459,187,134,36,2,91,11,324,18,3,9,341,8,71,67,17,37,396,59)), 
                  ncol=19,byrow=TRUE)
rownames(ctcf.pfm) <- c("A","C","G","T")

prior.params = c(A=0.25, C=0.25, G=0.25, T=0.25)
priorProbs = prior.params/sum(prior.params)
postProbs = t( t(ctcf.pfm + prior.params)/(colSums(ctcf.pfm)+sum(prior.params)) )
ctcf.pwm = Biostrings::unitScale(log2(postProbs/priorProbs))

library(BSgenome.Hsapiens.UCSC.hg19)
hg19 = BSgenome.Hsapiens.UCSC.hg19

p = patternMatrix(pattern=ctcf.pwm, windows=ctcf.peaks, genome=hg19, min.score=0.8)
heatMatrix(p, xcoords=c(-500, 500), main="CTCF motif")

plotMeta(mat=p, xcoords=c(-500, 500), smoothfun=function(x) stats::lowess(x, f = 1/10), 
         line.col="red", main="ctcf motif")