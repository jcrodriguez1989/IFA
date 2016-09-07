library(genefu);
library(limma);

##-----------------------------------------------------------------------------
##loadBCDataset high level constructor
##-----------------------------------------------------------------------------
##this function load the "genefu" data sets plus some other from Parket et al.
##It return a "limma" MAList object appropriate to be used by genefu package
loadBCDataset <- function(libname=c("upp", "nki", "vdx", "mainz", "transbig", "unt"), verbose=getOption("verbose", default=FALSE)){
    ##Check database
    stopifnot(!missing(libname))
    stopifnot(length(libname)==1)
    stopifnot(libname %in% c("upp", "nki", "vdx", "mainz", "transbig", "unt")) 

    ##Load the appropriate dataset and get its dataset name
    datasetName<-switch(libname,
	upp={require("breastCancerUPP"); "upp"},
	nki={require("breastCancerNKI"); "nki"},
	vdx={require("breastCancerVDX"); "vdx"},
	transbig={require("breastCancerTRANSBIG"); "transbig"},
	mainz={require("breastCancerMAINZ"); "mainz"},
	unt={require("breastCancerUNT"); "unt"})

    ##Data unification: annotation, clinical and expression.
    ##Load into memory and remove original variable
    if(verbose){message("Loading dataset...")}
    data(list=datasetName, envir=environment())
    dataset<-eval(parse(text=datasetName))
    rm(list=datasetName, envir=environment())

    ##Get annotation: probe, entrez and symbol
    if(verbose){message("Building the object")}
    annotData <-fData(dataset) 
    ##Check ProbeName field
    probe<-c("probe", "ProbeName")
    stopifnot(any(probe %in% names(annotData)))
    probe<-probe[probe %in% names(annotData)][1]
    ##Check EntrezID field
    entrez<-c("EntrezGene.ID", "EntrezID")
    stopifnot(any(entrez %in% names(annotData)))
    entrez<-entrez[entrez %in% names(annotData)][1]
    ##Check NCBI.gene.symbol field
    symbol<-c("Gene.symbol", "HUGO.gene.symbol")
    stopifnot(any(symbol %in% names(annotData)))
    symbol<-symbol[symbol %in% names(annotData)][1]

    ##Finally assign the fields
    annotData$probe <- as.character(annotData[,probe])
    annotData$EntrezGene.ID  <- as.integer(as.character(annotData[,entrez]))
    annotData$NCBI.gene.symbol <- as.character(annotData[, symbol])

    ##Get clinical and experimental data
    clinData <- pData(dataset)
    expData <- exprs(dataset) 

    ##Create PAM50 object
    MAaux<-new("MAList",list(M=expData,genes=annotData,targets=clinData))
    
    return(MAaux)
}

##-----------------------------------------------------------------------------
## PAM50Call
##-----------------------------------------------------------------------------
##Get PAM50 subtype 
##If std=="median" probes with the same mapping are averaged. This is done in 
##order to assure selecting the same "gene" to those in "genefu" library, 
##instead of the most variant probe (default in geneid.map), when more than one 
##probe match the same gene. This selection is based on probe population 
##variance that could depends on the number of accounted genes.
PAM50Call <- function(object, std=c("median", "none", "scale", "robust")[1], verbose=getOption("verbose", default=FALSE)){
    ##Check for standardization method
    stopifnot(std %in% c("median", "none", "scale", "robust"))
    stopifnot(length(std)==1)

    ##auxiliary pam
    pam50.aux <- pam50
    pam50.aux$std<-std

    ##dataset center scale using median estimation using all subjects present.
    if(std=="median"){
    ##Average probes with the same EntrezGene.ID
	if(verbose){message("Averaging over identical EntrezGeneID probes...")}
	object<-avereps(object,ID=as.character(object$genes$EntrezGene.ID))

	if(verbose){message("Center scaling using gene median...")}
	row.names(object$M) <- object$genes$probe
	med.total<-apply(object$M, 1, median, na.rm=TRUE)
	object$M.scale<-t(scale(t(object$M), center=med.total, scale=FALSE))
	row.names(object$M.scale) <- object$genes$probe
	dataset<-t(object$M.scale)
	pam50.aux$std<-"none"
    }else{
	dataset<-t(object$M)
    }

    ##Get the PAM50 calls
    if(verbose){message("Getting PAM50 subtypes...")}
    subtype.norm <- intrinsic.cluster.predict(sbt.model=pam50.aux, 
    data=dataset, annot=object$genes, do.mapping=TRUE, 
    do.prediction.strength=FALSE, verbose=verbose)

    ##Update object
    if(std != "median"){
	object$M.scale<-t(subtype.norm$profiles)
    }
    ##New fields
    object$subtype <- data.frame(subtype=factor(subtype.norm$subtype,
    levels=c("Basal", "Her2", "LumA", "LumB", "Normal")), 
    row.names=colnames(object$M))
    object$subtype.proba <-subtype.norm$subtype.proba 
    object$subtype.cor <- subtype.norm$cor 

    validObject(object)
    return(object)
}

