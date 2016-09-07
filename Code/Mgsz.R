require(mGSZ);

Mgsz <- function(exprMatrix, classes, genesets, GSEAcutoff) {
    mGSZ.obj <- mGSZ(x=exprMatrix, l=classes, y=genesets);
    mGszRes <- mGSZ.obj$mGSZ;
    
    mGszRes$Enriched <- mGszRes$pvalue <= GSEAcutoff;
    print(paste(sum(mGszRes$Enriched, na.rm=!F), " enriched terms"), sep="");
    return(mGszRes);
}
