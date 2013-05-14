makePxG <- function(cross, chr, marker, pheno.col) {
  pxg <- as.data.frame(cbind(geno=pull.geno(cross, chr=chr)[,marker], pheno=pull.pheno(cross, pheno.col=pheno.col)))
  pxg$geno <- factor(pxg$geno, labels=c("N2", "CB4856"))
  pxg$geno <- factor(pxg$geno, levels=c("N2", "CB4856"))
  #pxg <- droplevels(na.omit(pxg))
  return(pxg)
}

plotLOD <- function(scanone.obj, cutoff) {
  ggplot(data=as.data.frame(scanone.obj)) + aes(x=pos/1e6, y=lod) +
  geom_line(size=1) + geom_hline(data=as.data.frame(cutoff), linetype=2, aes(yintercept=cutoff)) +
  xlab("Genomic Position (Mb)") + ylab("LOD score") + facet_grid(.~chr, scales = "free_x", space="free") +
  scale_x_continuous(breaks=c(1:4*5)) + opts(legend.position="none")
}

plotLOD_talk <- function(scanone.obj, cutoff) {
  ggplot(data=as.data.frame(scanone.obj)) + aes(x=pos/1e6, y=lod) +
  geom_line(size=1) + geom_hline(data=as.data.frame(cutoff), linetype=2, aes(yintercept=cutoff)) +
  xlab("Genomic Position (Mb)") + ylab("LOD score") + facet_grid(.~chr, scales = "free_x", space="free") +
  scale_x_continuous(breaks=c(1:4*5)) + 
  opts(legend.position="none", axis.title.x = theme_text(vjust=0, size=24, face="bold"), 
         axis.title.y=theme_text(size=24, angle=90, vjust=0.25, face="bold"), 
       axis.text.x=theme_text(size=14, face="bold"), axis.text.y=theme_text(size=14, face="bold"))
}



plotLODbw <- function(scanone.obj, cutoff) {
  ggplot(data=as.data.frame(scanone.obj)) + aes(x=pos/1e6, y=lod) +
  geom_line(size=1) + geom_hline(data=as.data.frame(cutoff), linetype=2, aes(yintercept=cutoff)) +
  xlab("Genomic Position (Mb)") + ylab("LOD score") + facet_grid(.~chr, scales = "free_x", space="free") +
  scale_x_continuous(breaks=c(1:4*5)) + theme_bw() + opts(legend.position="none")
}

plotPxGbee <- function(pxg, marker, chr){
  ggplot(data=pxg) + aes(x=geno, y=pheno, color=geno, fill=geno) +
    geom_dotplot(binaxis="y", stackdir="center", method="histodot", binwidth=0.1) +
    xlab(paste("Genotype at marker", marker, "on chromosome", chr)) + ylab("RIAIL phenotype") +
    scale_color_manual(values=c("orange", "blue")) +
    scale_fill_manual(values=c("orange", "blue")) +
    opts(legend.position="none", axis.title.x = theme_text(vjust=0, size=18, face="bold"), 
       axis.title.y=theme_text(size=18, angle=90, vjust=0.25, face="bold"), 
       axis.text.x=theme_text(size=16, face="bold"), axis.text.y=theme_text(size=16, face="bold"))
}


plotPxGfull <- function(pxg, marker, chr) {
  ggplot(data=pxg) + geom_jitter(size=2, position=position_jitter(width=0.15), aes(x=geno, y=pheno, colour=geno)) +
  geom_boxplot(alpha=0.5, width=0.3, outlier.size=0, aes(x=geno, y=pheno, fill=geno)) + 
  xlab(paste("Genotype at marker", marker, "on chromosome", chr)) + ylab("RIAIL phenotype") +
  scale_color_manual(values=c("orange", "blue")) +
  scale_fill_manual(values=c("orange", "blue")) + opts(legend.position="none")
}

plotPxGbox <- function(pxg, marker, chr) {
  ggplot(data=pxg) + geom_boxplot(alpha=0.5, width=0.3, outlier.size=0, aes(x=geno, y=pheno, fill=geno)) + 
  xlab(paste("Genotype at marker", marker, "on chromosome", chr)) + ylab("RIAIL phenotype") +
  scale_color_manual(values=c("orange", "blue")) +
  scale_fill_manual(values=c("orange", "blue")) + opts(legend.position="none")
}

plotPxGdots <- function(pxg, marker, chr) {
  ggplot(data=pxg) + geom_jitter(size=4, position=position_jitter(width=0.15), aes(x=geno, y=pheno, colour=geno)) +
  xlab(paste("Genotype at marker", marker, "\non chromosome", chr)) + ylab("RIAIL phenotype") +
  scale_color_manual(values=c("orange", "blue")) +
  scale_fill_manual(values=c("orange", "blue")) + 
  opts(legend.position="none", axis.title.x = theme_text(vjust=0, size=24, face="bold"), 
       axis.title.y=theme_text(size=24, angle=90, vjust=0.25, face="bold"), 
       axis.text.x=theme_text(size=14, face="bold"), axis.text.y=theme_text(size=14, face="bold"))
}

  #Need to fix
plotDensity <- function(pxg) {
  ggplot(data=pxg) + geom_density(size=2, aes(x=pheno, colour=geno)) +
  xlab("RIAIL phenotype") + ylab("Statistical density")
  scale_color_manual(values=c("orange", "blue")) +
  opts(legend.position="none")
}
  
 #Need to fix
plotPhenoDist <- function(pxg, phenotype) {
  ggplot(data=pxg) + aes(x=pheno) + geom_histogram() + 
  xlab(paste(phenotype)) + ylab("Count") + 
  opts(legend.position="none", axis.title.x = theme_text(vjust=0, size=24, face="bold"), 
       axis.title.y=theme_text(size=24, angle=90, vjust=0.25, face="bold"), 
       axis.text.x=theme_text(size=14, face="bold"), axis.text.y=theme_text(size=14, face="bold"))
  
}

plotPxGparents <- function(pxg, parents, pheno.col, phenotype) {
  pxg$riail <- 1
  temp.riails <- data.frame(pheno = pxg$pheno, riail= pxg$riail)
  temp.n2 <- data.frame(pheno = parents[parents$strain=="N2", pheno.col], riail = parents$riail[parents$strain=="N2"])
  temp.cb <- data.frame(pheno = parents[parents$strain=="CB4856", pheno.col], riail = parents$riail[parents$strain=="CB4856"])
  temp <- rbind(temp.n2, temp.cb, temp.riails)
  temp$riail <- factor(temp$riail, levels=c(3, 4, 1), labels=c("N2", "CB4856", "RIAILs"))
  ggplot(data=temp) + aes(x=riail, y=pheno, colour=riail) + geom_jitter(size=2, position=position_jitter(width=0.15)) +
    xlab("") + ylab(paste(phenotype)) + scale_color_manual(values=c("orange","blue","black")) + opts(legend.position="none")
}

#Need to normalize for different density peaks

plotPxGparents_dense <- function(pxg, parents, pheno.col, phenotype) {
  pxg$riail <- 1
  temp.riails <- data.frame(pheno = pxg$pheno, riail= pxg$riail)
  temp.n2 <- data.frame(pheno = parents[parents$strain=="N2", pheno.col], riail = parents$riail[parents$strain=="N2"])
  temp.cb <- data.frame(pheno = parents[parents$strain=="CB4856", pheno.col], riail = parents$riail[parents$strain=="CB4856"])
  temp <- rbind(temp.n2, temp.cb, temp.riails)
  temp$riail <- factor(temp$riail, levels=c(3, 4, 1), labels=c("N2", "CB4856", "RIAILs"))
  ggplot(data=temp) + aes(x=pheno, colour=riail) + geom_density() +
    xlab(paste(phenotype)) + ylab("Statistical density") + scale_color_manual(values=c("orange","blue","black")) + opts(legend.position="none")
}

