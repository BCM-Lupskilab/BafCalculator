getAOH_data <- function(snpDir,sampleName)
{
	
	source ("vcf_utils.R")
	fileSNP <- paste(snpDir,"/",sampleName,".SNPs_Annotated.vcf.bz2", sep="")
	dataSNP <- read.vcf.quick(fileSNP)
	
	# parsing genotype
	gt<- do.call(rbind, strsplit(dataSNP[[sampleName]], ":"))
	vR <- gt[,2]
	tR <- gt[,4]
	dataSNP[,vR:=as.numeric(vR)]
	dataSNP[,tR:=as.numeric(tR)]
	dataSNP[,af:=(as.numeric(vR)/ as.numeric(tR))]
	
	
	# segmentation
	dataSNP <- dataSNP[which(dataSNP$FILTER == "PASS"),]
	dataSNP[,baf:=abs(af - 0.5)]
	library("DNAcopy")
	CNA.obj = CNA(dataSNP$baf, dataSNP$CHROM, dataSNP$POS, data.type = "binary")

	segment.obj = segment(CNA.obj,   verbose = 1)
	out <- segment.obj$output
	out
    
	xx <- dataSNP
	return(list(xx =xx, out=out))
}




plotSample<- function(sampleName, snpDir,  chrom, start,stop,stopChrom, wd, mutation, allExons, allGenes, all=F, threshold=0.47 ,sizeThreshold=0.5e6)
{
    chr  <- paste("chr", chrom, sep="")
    wdR <- wd
    wdL <- wd
    
    # reading snp data
    fileSNP <- paste( snpDir,'/',sampleName,".SNPs_Annotated.vcf.bz2", sep="")
    dataSNP <- read.vcf.quick.noinfo(fileSNP)
    
    # parsing genotype
    gt<- do.call(rbind, strsplit(dataSNP[[sampleName]], ":"))
    vR <- gt[,2]
    tR <- gt[,4]
    dataSNP[,vR:=as.numeric(vR)]
    dataSNP[,tR:=as.numeric(tR)]
    dataSNP[,af:=(as.numeric(vR)/ as.numeric(tR))]
    
    
    # segmentation
    dataSNP <- dataSNP[which(dataSNP$FILTER == "PASS"),]
    dataSNP[,baf:=abs(af - 0.5)]
    library("DNAcopy")
    CNA.obj = CNA(dataSNP$baf, dataSNP$CHROM, dataSNP$POS, data.type = "binary")
    segment.obj = segment(CNA.obj,   verbose = 1)
 
    out <- segment.obj$output
    xx <- dataSNP
    
    
    
    plotGenes<- FALSE
    cexv<- 1.5
    
    if (all == TRUE)
    {
        start <- 0
        stop <- stopChrom
        plotGenes <- FALSE
        cexv<- 0.5
        wd <- 0
        wdR <- 0
        wdL <- 0
    }
    
    ex <- allExons[which((allExons$Chr == chr) &(allExons$Stop >= (start -wdL)) &( allExons$Start < (stop + wdR))),]
    gn <- allGenes[which((allGenes$Chr ==chr) &(allGenes$Stop >= (start -wdL))& (allGenes$Start < (stop + wdR))),]
    if (length(which(duplicated(gn$Name))))
    {
        gn <- gn[-which(duplicated(gn$Name)),]
    }

    
    library(Hmisc)
    dd <- xx[which((xx$CHROM==chrom) ),]
    dd <- dd[order(dd$CHROM, as.numeric(dd$POS)),]
    
    minY <- -0.4
    if (!plotGenes){minY <- 0}
    
    ## display empty plot
    plot(dd$POS, dd$af, ylim=c(minY, 1.1), xaxs = "i", cex=cexv,yaxt="n", xlab="", ylab="variants/total reads ratio", xlim=c(start-wdL,stop+wdR), type="n",cex.lab=1.5,cex.axis=1.5)
    
    
    ## plot predicted segments
    segData <- out[which(out[,2]==chrom),]
    selData <- segData[which(segData[,6] > threshold & (segData[,4] - segData[,3]) > sizeThreshold),]
    if (!all){print(selData)}
    lapply(1:nrow(selData), function(i){
        row <- selData[i,]
        rect(row[3],0, row[4], 1, col="lightgray" , border=NA)
    })
    
    
    ## plot SNP data
    points(dd$POS, dd$af,pch=16, cex=cexv,yaxt="n")
    axis(2, at=1:10/10,labels=as.character(1:10/10), col.axis="black", las=2)
    abline(v=c(mutation), col="red",lwd=2)
    #abline(v=c(mutation1), col="blue",lwd=2)
    
    ## plot genes
    if (plotGenes)
    {
        par(lend=2);aa <- lapply(1:nrow(gn), function(i){
            offset <- 0.1 +i%%5 * 0.05
            g <-gn[i,] ;    lines(c(g$Start, g$Stop),c(0-offset,0-offset), col="blue")
            ee <- ex[ex$Name == g$Name,]
            lapply(1:nrow(ee), function(i){e <-ee[i,] ;    rect(e$Start,0.02 - offset, e$Stop, -0.02-offset, col="darkblue", border="darkblue")})
        })
        
        abline(h= 0,col="lightgray", lwd=1)
        
    }
    
    
}
	addColNames2lff <- function (x)
	{
		colnames(x) <- c("ClassName", "Name", "Type", "Subtype", "Chr", "Start" , "Stop" , "Strand", "Phase", "Score" , "Other" )
		return (x[,c("ClassName", "Name", "Type", "Subtype", "Chr", "Start" , "Stop" , "Strand", "Phase", "Score" , "Other" )])
	}
	
	
############### USAGE ###########
	
	

allExons <- addColNames2lff(read.csv("allExons.lff", stringsAsFactors=F, header=F, sep="\t"))
allGenes <- addColNames2lff(read.csv("allGenes.lff", stringsAsFactors=F, header=F, sep="\t"))

sampleName <- 'BAB5133' # sampleName

snpDir <- paste("/home/va/data/raw_data/",sampleName,sep="") # directory with VCF file
threshold <- 0.47 # threshold for segmentation
sizeThreshold <- 0.5e6 # threshold for segmentation


chrom <- "9" # chromsome to display
start <-  116918231 # start gene pos
stop<- 117072975 # stop gene pos
stopChrom <-  141213431  # stop pos of chromosome
wd <- 24e6 # left and right extension

mutation <- 116999951   # pos of mutation


outputFilename <- paste("~/BAB5133-chr9.pdf", sep="") # output file names
pdf (file=outputFilename, width=1.5*10, height=0.75*12.5) ## Plotting SVG

#### plotting
plotSample(sampleName, snpDir,  chrom, start,stop,stopChrom, wd, mutation,allExons, allGenes, all=T, threshold=threshold, sizeThreshold=sizeThreshold)
plotSample(sampleName, snpDir,  chrom, start,stop,stopChrom, wd, mutation,allExons, allGenes, all=F, threshold=threshold, sizeThreshold=sizeThreshold)
dev.off()





