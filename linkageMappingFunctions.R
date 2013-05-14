library(qtl)
#getting the cross object with the ready-made geno part, the mapping and a couple useful functions
load("~/Dropbox/Coding/Active Projects/Report 2/N2xCB4856_RIAILs_Rqtlfiles.RData")
genomic.pos <- read.csv(file="~/Dropbox/Coding/Active Projects/Report 2/SNPgenomicpositions.csv", header=TRUE, stringsAsFactors=FALSE)
names(genomic.pos) <- c("gPos", "marker.name")

# Just to make things a little neater; we can turn this off if there are problems later
options(warn=-1)

######### The following code was taken and adapted from Josh Bloom (full code can be found in Active Projects/ecaPeakCode.R)  ###############

##### Extract genotype matrix from cross structure and recode as -1,1 #####################################
    # originally stored as 1s and 2s in cross object
extractScaledGenotype=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$data }))*2)-3 }

#alternatively, don't change the 1s and 2s in the cross object. 
#NOTE: this means you want intercept to be TRUE in getPhenoResids if this is the argument you provide as gdata
extractGenotype=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$data })))}
   
#takes in the phenotype data.frame from Report 1 and scales the data if desired; have to exclude id column before
#combining back together. Optional scaleVar argument also taken in
scalePhenotype=function(pheno.data.frame,scaleVar=FALSE){cbind("id"=pheno.data.frame$id, apply(pheno.data.frame[colnames(pheno.data.frame)!="id"], 2, scale, scale=scaleVar))}

#gets the maximum LOD score of the scanone object for each chromosome
#also gets the marker index (not the position) that will later be matched up with gdata in getPhenoResids
getChrPeaks = function(yourScanOne) {
		###getting the max LOD values per chromosome
	    chr.peaks.lod = tapply(yourScanOne$lod, yourScanOne$chr, max)
	    ###finding the number of markers on each chromosome
        chr.mindex.offset = tapply(yourScanOne$lod, yourScanOne$chr, length)
        ###uses the number of markers per chromosome to find the offset for the max LOD marker index.
        ###when finding the max LOD marker index, which.max will start at 1 for each chromosome,
        ###but we need it to match up to gdata later, and the offset will allow that.
        chr.mindex.offset = c(0, cumsum(chr.mindex.offset)[1:5])
        ### get marker index of LOD peaks per chromosomes                             
        chr.peaks.index = tapply(yourScanOne$lod, yourScanOne$chr, which.max)
        ### add offset to LOD peak marker index                             
        chr.peaks.index = chr.peaks.index + chr.mindex.offset
		#combines the peak LOD numbers and their marker indices
        return(list(chr.peaks.lod = chr.peaks.lod, chr.peaks.index=chr.peaks.index))
    }
    
#filters the LOD peaks from getChrPeaks based on the threshold you give it
getPeakArray = function(peaklist, threshold, multiTraits=FALSE) {
        tryCatch( {
        keepPeaks   = which(peaklist$chr.peaks.lod>threshold, arr.ind=T)
        if(multiTraits){
		        kP = data.frame(rownames(keepPeaks), peaklist$chr.peaks.index[keepPeaks])
   		        names(kP)=c('trait', 'markerIndex') 
        		kP = kP[order(kP$trait, kP$markerIndex),]
        		return(kP)} 
        else{
        		kP = data.frame(peaklist$chr.peaks.index[keepPeaks])
        		names(kP)="markerIndex"
        		return(kP)}
        		}, error=function(e) {return(NULL) })		
    }
    
#regresses out the peaks given by getPeakArray
#cross.phenotypes is your cross object$pheno, gdata is either scaled or unscaled
#IMPORTANT: if gdata is unscaled (data is recorded as 1s and 2s), intercept=TRUE
getPhenoResids = function(cross.phenotypes,gdata, peakArray, intercept=FALSE) {
        presids = cross.phenotypes
        spA = peakArray
        if(length(spA)>0){
            if(intercept) {
                 rr = residuals(lm(cross.phenotypes[,yourTransformation]~gdata	[,spA$markerIndex]))
            }else{
                 rr = residuals(lm(cross.phenotypes[,yourTransformation]~gdata[,spA$markerIndex]-1))
            }
            presids[as.numeric(names(rr)),yourTransformation]=rr
        }
        return(presids)
    }
    
    
 doQTLModel = function(peak.index, fcross, LODSso, refine=TRUE){
 	### As before, yourTransformation defined in R Markdown file by user
        qtlMarkers = do.call('rbind', peak.index)$markerIndex
        qtlMarkers = unique(qtlMarkers)
        mqtl=NA
        fqtl=NA
        aqtl=NA
        rqtl=NA
        qtlMarkersCP = NA
        CIs =NA
        nQTL = length(qtlMarkers)
        if(nQTL>0) {
            qtlMarkersNames = rownames(LODSso)[qtlMarkers]
            qtlMarkersCP = LODSso[qtlMarkersNames ,c(1,2)]

            mqtl =  makeqtl(fcross, chr=qtlMarkersCP$chr, pos=qtlMarkersCP$pos, qtl.name=qtlMarkersNames, what='prob')
            if(refine == TRUE) {
                rqtl  = refineqtl(fcross,pheno.col=yourTransformation, mqtl, method='hk',
                               model='normal', incl.markers=TRUE,maxit=1, verbose=FALSE)
            } else {
                 rqtl = mqtl
            }
            fqtl  = fitqtl(fcross, pheno.col =yourTransformation, rqtl, method='hk', get.ests=TRUE, dropone=TRUE)
                
            #scan for interactors amongst QTLs with significant main effects
             if(nQTL==1) { aqtl=NA } 
             if(nQTL>1)  { aqtl = addint(fcross, pheno.col=yourTransformation, rqtl, method='hk', qtl.only=TRUE)      }

            # get confidence intervals
            if(refine ==TRUE) {
                CIs = lapply(1:length(qtlMarkers),function(x) {lodint(rqtl,qtl.index=x) })
                }
        }
       # we want to keep the scanone object so it  an feed into makeQTL.df.keep.so
       return(list(qtlMarkers=qtlMarkers, qtlMarkersCP=qtlMarkersCP, mqtl=mqtl, fqtl=fqtl, aqtl=aqtl,rqtl=rqtl,CIs=CIs))
    } 


#putting all the peak finding parts together now...
findAllThePeaks = function(impcross, phenotypes, scaleGeno=FALSE, scalePheno=TRUE, model=c("normal", "np")) {
		if(model=="np"){yourTranformation <- 5}
		impcross <- calc.genoprob(impcross, step=0)
		if(scaleGeno){
			gdata <- extractScaledGenotype(impcross)
		}else{
			gdata <- extractGenotype(impcross)
		}
		if(scalePheno){
			pdata <- scalePhenotype(phenotypes)
		}else{
			pdata <- phenotypes
		}
		#genotype marker framework already in place, just need to add in the phenotypes
		#mergePheno defined from loading N2xCB4856_RIAILs_Rqtlfiles.RData
		impcross$pheno <- mergePheno(impcross, pdata)
		#yourTransformation defined in the R Markdown file
		mapped.trait <- scanone(impcross, pheno.col=yourTransformation, model=model, method="mr", verbose=FALSE)
		#finds the LOD score significance threshold based on 1000 permutations of the data
		threshold.mapping <- scanone(impcross, pheno.col=yourTransformation, model=model, method="mr", n.perm=1000, verbose=FALSE)
		#cutoff for LOD significace is the 95th quantile of the permutation mapping above
		threshold.cutoff <- quantile(threshold.mapping, probs=0.95)[[1]]
		#getting max LOD score per chromosome
		LOD.peaks <- getChrPeaks(mapped.trait)
		#filters out LOD.peaks and keeps only the significant ones
		LOD.peaks.significant <- getPeakArray(LOD.peaks, threshold.cutoff)
		#until we've found all the significant LOD peaks, this train will keep chuggin'
		iterations = 1
		allThePeaks <- list(LOD.peaks.significant)
		regressed.cross <- impcross
		while(TRUE){
			iterations <- iterations + 1
			if(scaleGeno){
				regressed.cross$pheno <- getPhenoResids(regressed.cross$pheno, gdata, LOD.peaks.significant, intercept=FALSE)
				}
			else{
				regressed.cross$pheno <- getPhenoResids(regressed.cross$pheno, gdata, LOD.peaks.significant, intercept=TRUE)
				}
			new.mapped.trait <- scanone(regressed.cross, pheno.col=yourTransformation, model=model, method="mr")
			new.LOD.peaks <- getChrPeaks(new.mapped.trait)
			LOD.peaks.significant <- getPeakArray(new.LOD.peaks, threshold.cutoff)
			if(nrow(LOD.peaks.significant)==0){break}
			else{allThePeaks[iterations] <- list(LOD.peaks.significant)}
		}
		allQtl <- doQTLModel(allThePeaks, impcross, mapped.trait, refine=TRUE)
		# we'll need all of the information here for to make plots and get the dataframe
		return(list(QTLs=allQtl, cutoff=threshold.cutoff, SO=mapped.trait))
}

makeQTL.df = function(QTL.list, ScanOne.object, mappingName){
	totalEffectSize <- round(QTL.list$fqtl$result.full["Model", "%var"], digits=3)
	if(length(QTL.list$qtlMarkers)>1){
			percentVar <- round(QTL.list$fqtl$result.drop[, "%var"], digits=3)
	}
	else{
			percentVar <- totalEffectSize
		}
	leftCPos <- sapply(1:length(QTL.list$CIs), function(x){QTL.list$CIs[[x]]$pos[1]})
	rightCPos <- sapply(1:length(QTL.list$CIs), function(x){QTL.list$CIs[[x]]$pos[3]})
	leftGPos.marker.names <- sapply(1:length(QTL.list$CIs), function(x){rownames(QTL.list$CIs[[x]])[1]})
	leftGPos <- genomic.pos[match(leftGPos.marker.names, genomic.pos$marker.name), 1]
	rightGPos.marker.names <- sapply(1:length(QTL.list$CIs), function(x){rownames(QTL.list$CIs[[x]])[3]})
	rightGPos <- genomic.pos[match(rightGPos.marker.names, genomic.pos$marker.name), 1]
	peakGPos.marker.names<- sapply(1:length(QTL.list$CIs), function(x){rownames(QTL.list$CIs[[x]])[2]})
	peakGPos <- genomic.pos[match(peakGPos.marker.names, genomic.pos$marker.name), 1]
	peakCPos <- ScanOne.object$pos[QTL.list$qtlMarkers]
	LODscore <- ScanOne.object$lod[QTL.list$qtlMarkers]
	chr <- ScanOne.object$chr[QTL.list$qtlMarkers]
	QTL.df <- data.frame("mapping"=mappingName, "VE"=percentVar, "chr"=chr, "leftCPos"=leftCPos, "rightCPos"=rightCPos, "peakCPos"=peakCPos, "LODscore"=LODscore, "leftGPos"=leftGPos, "rightGPos"=rightGPos, "peakGPos"=peakGPos)
	return(QTL.df)
}
    
plotQTLs = function(QTL.df, ScanOne.object, threshold, model){
	#genomic.pos defined at the top of this file, at lines 4 and 5
	ScanOne.object <- data.frame(ScanOne.object, "gPos"=genomic.pos$gPos)
	plot <- ggplot(ScanOne.object) + geom_line(size=1) + aes(x=gPos/1e+06, y=lod) + 
	geom_hline(yintercept=threshold, linetype="dashed") + 
    geom_point(data=QTL.df, aes(x=peakGPos/1e+06, y=LODscore + 0.3), colour="red", shape=8, size=5) +
    geom_segment(data=QTL.df, aes(xend=rightGPos/1e+06, x=leftGPos/1e+06, y=LODscore + 0.1, yend=LODscore+0.1), colour="red", linetype="dashed") +
    geom_text(data=QTL.df, aes(x=peakGPos/1e+06, y=LODscore+0.5, label=VE), size=4) +
    labs(title=paste("Peaks Found With", model, "Mapping"), x="Genomic Position (Mb)", y="LOD score") +
    facet_grid(.~chr, scales="free_x") +
    theme(title=element_text(face='bold', size=20), axis.text=element_text(face='bold', size=12))
}

makePxG <- function(impcross, pheno.data.frame, markerIndex, scalePheno=TRUE) {
  if(scalePheno){
			pdata <- scalePhenotype(pheno.data.frame)
  }else{
			pdata <- pheno.data.frame
  }
  impcross$pheno <- mergePheno(impcross, pdata)
  gdata <- extractGenotype(impcross)
  QTLgdata <- gdata[,markerIndex]
  QTLpdata <- impcross$pheno[, yourTransformation]
  pxg <- data.frame("geno"=QTLgdata, "pheno"=QTLpdata)
  pxg$geno <- factor(pxg$geno, labels=c("N2", "CB4856"))
  pxg$geno <- factor(pxg$geno, levels=c("N2", "CB4856"))
  #pxg <- droplevels(na.omit(pxg))
  return(pxg)
}

plotPxGfull <- function(pxg, QTL.df, QTL.index, model) {
  marker = rownames(QTL.df)[QTL.index]
  chr = QTL.df$chr[QTL.index]
  ggplot(data=pxg) + geom_jitter(size=2, position=position_jitter(width=0.15), aes(x=geno, y=pheno, colour=geno)) +
  geom_boxplot(alpha=0.5, width=0.3, outlier.size=0, aes(x=geno, y=pheno, fill=geno)) + 
  xlab(paste("Genotype at marker", marker, "on chromosome", chr)) + ylab("RIAIL phenotype") +
  scale_color_manual(values=c("orange", "blue")) +
  scale_fill_manual(values=c("orange", "blue")) + theme(legend.position="none") +
  labs(title=paste("Phenotype at QTL Split by Strain Using", model, "Mapping"))
}

merged.cross <- N2xCB4856.cross
merged.cross$pheno <- mergePheno(N2xCB4856.cross, yourPhenoTransforms.df)

normalQTL.data <- findAllThePeaks(N2xCB4856.cross, yourPhenoTransforms.df, model="normal")
normalQTL.list <- normalQTL.data$QTLs
normalQTL.cutoff <- normalQTL.data$cutoff
normalQTL.SO <- normalQTL.data$SO
normalQTL.model <- makeQTL.df(normalQTL.list, normalQTL.SO, "normal")
normalQTL.plot <- plotQTLs(normalQTL.model, normalQTL.SO, normalQTL.cutoff, "Normal")
normalQTL.RIAIL.distro <- lapply(1:nrow(normalQTL.model), function(x){plotPxGfull(makePxG(N2xCB4856.cross, yourPhenoTransforms.df, normalQTL.list$qtlMarkers[x]), normalQTL.model, x, "Normal")})
normalQTL.effects <- effectscan(sim.geno(merged.cross, step=0, stepwidth="variable"), pheno.col=5, draw=FALSE)
normalQTL.effects <- ggplot(normalQTL.effects) + geom_line(size=1) + aes(x = pos, y = a) + facet_grid(.~chr, scales="free_x") + geom_vline(data=normalQTL.model, aes(xintercept=peakCPos), colour="red", linetype="dashed") + labs(x = "Genetic Position", y = "% VE Explained", title="Variant Effect Scan") + theme(title=element_text(face='bold', size=20), axis.text=element_text(face='bold', size=12))
#print(normalQTL.model)
#print(normalQTL.plot)
#print(normalQTL.RIAIL.distro)

npQTL.data <- findAllThePeaks(N2xCB4856.cross, yourPhenoTransforms.df, model="np")
npQTL.list <- npQTL.data$QTLs
npQTL.cutoff <- npQTL.data$cutoff
npQTL.SO <- npQTL.data$SO
npQTL.model <- makeQTL.df(npQTL.list, npQTL.SO, "nonparametric")
npQTL.plot <- plotQTLs(npQTL.model, npQTL.SO, npQTL.cutoff, "Nonparametric")
npQTL.RIAIL.distro <- lapply(1:nrow(npQTL.model), function(x){plotPxGfull(makePxG(N2xCB4856.cross, yourPhenoTransforms.df, npQTL.list$qtlMarkers[x]), npQTL.model, x, "Nonparametric")})
#print(npQTL.model)
#print(npQTL.plot)
#print(npQTL.RIAIL.distro)

cim.hk <- cim(merged.cross, pheno.col=yourTransformation, method="hk")
cim.plot <- ggplot(cim.hk) + geom_line() + aes(x=pos, y=lod) + facet_grid(.~chr, scales="free_x") + labs(title="Composite Interval Mapping", x="Genetic Position", y="LOD score")
if(cimPermutations){
		cim.cutoff <- cim(merged.cross, pheno.col=yourTransformation, method="hk", n.perm=1000)
		cim.cutoff <- quantile(cim.cutoff, probs=0.95)[[1]]
		cim.plot <- ggplot(cim.hk) + geom_line() + aes(x=pos, y=lod) + facet_grid(.~chr, scales="free_x") + labs(title="Composite Interval Mapping", x="Genetic Position", y="LOD score") + geom_hline(yintercept=cim.cutoff, linetype="dashed")
		}


#print(cim.plot)

overlaid.plots <- ggplot() + geom_line(data=normalQTL.SO, colour="red") + aes(x=pos, y=lod) + facet_grid(.~chr) + geom_line(data=npQTL.SO, colour="brown") + 
					geom_line(data=cim.hk, colour="blue") + labs(title="Overlaid Plots for Normal, Nonparametric and Composite Interval Mapping", x="Genetic Position", y="LOD score") +
					geom_point(data=normalQTL.model, aes(x=peakCPos, y=LODscore + 0.4), colour="red", shape=8, size=5) + geom_point(data=npQTL.model, aes(x=peakCPos, y=LODscore + 0.4), colour="brown", shape=8, size=5)
#print(overlaid.plots)



	

