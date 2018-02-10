	
 getKVdf <- function(xx, split.char=",",removequotes=F)
{
	yy <- strsplit(xx, split.char)[[1]]
	kv <- do.call(rbind, lapply ( (strsplit(yy, "=")), function(ll){cbind(ll[1], paste(ll[2:length(ll)], collapse="="))}))
	colnames(kv) <- c("key", "value")
	if (removequotes){
		kv <- gsub ( "\"", "", kv)
	}
	kv <- data.frame(kv, stringsAsFactors=F)
	colnames(kv) <- c("key", "value")
	kv
}
read.vcf.header <- function(file){
		
		line <- scan(file=pipe(paste ("bzcat", file)), n=1, skip=0, sep="\n", what="character", quiet=TRUE)
		skip <- 1
		header <- c()
		while(substr(line, 0 ,1) == "#" ){
			header <- c(header, gsub("#", "", line) )
			line <- scan(file=file, n=1, skip=skip, sep="\n", what="character", quiet=TRUE)
			skip <- skip + 1
		}
		lines <- header
		
		getSubHeader <- function(headerTitle, lines){
			starts <- regexpr(headerTitle, lines) + nchar(headerTitle) 
			stops <- regexpr( ">", lines) -1
			
			hd <- lapply(1:length(starts), function(i){
						#			print (i)
						start <- starts[i]
						stop <- stops[i]
						if (start > nchar(headerTitle) )
						{
							xx <- substr(lines[i], start, stop)
							kv <- getKVdf(xx, removequotes=T)
						}
						
						else {NULL}
					})
			
			hdf <- Filter(function(x) !is.null(x), hd)
			names(hdf) <- lapply(hdf, function(x)x$value[which(x$key=="ID")])
			hdf
		}
		info <- getSubHeader( "INFO=<", lines)
		format <-  getSubHeader( "FORMAT=<", lines)
		
		list (header=strsplit(lines[length(lines)], "\t")[[1]], INFO=info,FORMAT =format, nlines=length(lines) + 1)
	}

	splitColNames <- function(data, col, header){
			xx <- lapply(names(	header[[col]]), function(nm){
						print (nm)
						pos <- regexpr(paste(nm,"=[^;]*;", sep=""), data[[col]])
						starts <- pos + nchar(nm) +1
						stops <- pos + attr(pos, "match.length") -2
						res <- substr(data[[col]],starts,  stops)
						res [which(starts <  0)] <- NA
						print (length(which(starts <  0)))
						data[, eval(nm) :=res]			
					})
		}
		
		splitColNames2 <- function(data, commonHeader, file){
			
			
			write (data$INFO, paste(file, ".INFO", sep=""))
			info <- fread (paste(file, ".INFO", sep=""), sep=";", header=FALSE)
			xx <- lapply(1:length(commonHeader), function(i){
						nm <- commonHeader[i]
						res <- substr(info[[i]], nchar(nm) + 2, 100)
						data[, eval(nm) :=res]
					}
			)
			system(paste ("rm " ,file, ".INFO", sep="") )
		}
		read.vcf <- function(file){
			header <- read.vcf.header(file)
			library(data.table)
			system (paste("bzcat ",file," | tail -n +",header$nlines," > " , file, ".noheader",  sep="" ))
			data <- fread(paste(file , ".noheader", sep=""), header=FALSE, stringsAsFactors=F, sep="\t")
			system(paste ("rm " ,file, ".noheader", sep="") )
			setnames(data, header$header)	
			commonHeader <- sapply ( strsplit(data$INFO[1], ";")[[1]], function(x){strsplit(x, "=")[[1]][1]})
			#splitColNames2(data, commonHeader, file)
			splitColNames(data, "INFO", header)
			setnames(data, commonHeader, sapply (header$INFO[commonHeader], function(x){gsub(" ","_", x$value[which(x$key =="Description")])}))
			data
		}
		
		
		read.vcf.quick <- function(file){
			header <- read.vcf.header(file)
			library(data.table)
			system (paste("bzcat ",file," | tail -n +",header$nlines," > " , file, ".noheader",  sep="" ))
			data <- fread(paste(file , ".noheader", sep=""), header=FALSE, stringsAsFactors=F, sep="\t")
			system(paste ("rm " ,file, ".noheader", sep="") )
			setnames(data, header$header)	
			#commonHeader <- sapply ( strsplit(data$INFO[1], ";")[[1]], function(x){strsplit(x, "=")[[1]][1]})
			#splitColNames2(data, commonHeader, file)
			#splitColNames2(data, "INFO", header)
			#setnames(data, commonHeader, sapply (header$INFO[commonHeader], function(x){gsub(" ","_", x$value[which(x$key =="Description")])}))
			data
		}
		
		read.vcf.quick.noinfo <- function(file){
			header <- read.vcf.header(file)
			library(data.table)
			system (paste("bzcat ",file," | tail -n +",header$nlines," | cut -f1-7,9,10 > " , file, ".noheader",  sep="" ))
			
			data <- fread(paste(file , ".noheader", sep=""), header=FALSE, stringsAsFactors=F, sep="\t")
			system(paste ("rm " ,file, ".noheader", sep="") )
			setnames(data, header$header[-9])	
			data
		}
		
		#fileSNP <- "/home/tgambin/workspace/mendelian-repo/tgambin/data/raw_data/11-0015/11-0015.SNPs_Annotated.vcf"
#fileSNP <- "/home/tgambin/workspace/mendelian-repo/tgambin/data/raw_data/BAB3986/BAB3986.SNPs_Annotated.vcf"
#fileINDEL <- "/home/tgambin/workspace/mendelian-repo/tgambin/data/raw_data/BAB3986/BAB3986.INDELs_Annotated.vcf"
#sampleName <- "BAB3986"

#fileSNP <- "/home/tgambin/workspace/mendelian-repo/tgambin/data/raw_data/548810/548810.SNPs_Annotated.vcf"
#fileINDEL <- "/home/tgambin/workspace/mendelian-repo/tgambin/data/raw_data/548810/548810.INDELs_Annotated.vcf"
#sampleName <-"548810"
#mergeAndConvertVCFs(fileSNP, fileINDEL, sampleName)

		mergeAndConvertVCFs <- function(snp, indel, sampleName)
		{
			#fileName <- paste(path , sampleName,sep="")
			fileINDEL <- indel
			fileSNP <- snp
			dataINDEL <- read.vcf(fileINDEL)
			dataSNP <- read.vcf(fileSNP)
			zz <- lapply (setdiff(names(dataINDEL),names(dataSNP)), function(nm){dataSNP [, eval(nm):="."]})
			zz <- lapply (setdiff(names(dataSNP),names(dataINDEL)), function(nm){dataINDEL [, eval(nm):= "."]})
			
			data <- rbind(dataSNP, dataINDEL[,names(dataSNP), with=F])
			toRM <- which(colnames(data) %in% c("INFO", "ID" ))
			dataTmp <- data[,-toRM,with=F]
			
			
			# put gene name to 11-th place - required by other after burners
			geneNamePos <- which(colnames(dataTmp) == "Gene name")
			dataFinal <- dataTmp[,c(1:10, geneNamePos , setdiff(11:(ncol(dataTmp)), geneNamePos)),with=F]
			
			tab<- strsplit (snp , "/")[[1]]
			path <- paste (tab[1:(length(tab) - 1)], collapse="/")
			
			write.table(dataFinal,  row.names=F,  file=paste(path,"/", sampleName, ".tsv", sep=""), sep="\t" , quote=F)
		}
		
		