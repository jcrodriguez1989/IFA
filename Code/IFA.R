source("DEnricher.R");
source("loadGO.R");
source("Mgsz.R");

## IFA executes Functional Analysis through mGSZ and dEnricher. If correctly set
## dEnricher parameters in ... then it will execute the desired dEnricher 
## analysis, and mGSZ with the same gene sets.
## exprMatrix: expression matrix. Genes as rows, samples as columns.
## Rownames must be EntrezGene IDs.
## classes: string vector of classes representing each column from exprMatrix.
## No more than two classes are allowed, classes length must be the same as
## exprMatrix's number of columns.
## genesets (optional): list of gene sets (each gene set is a vector of strings)
## If dEnricher parameters are correctly set in ... then dEnricher gene sets
## will be used.
## If not set, then the last Gene Ontology is loaded.
## SEAcutoff: enrichment cutoff value for SEA.
## GSEAcutoff: enrichment cutoff value for GSEA.
## br: choosen br for SEA analysis, it can be a vector of Entrez IDs to be used
## as background reference. It also can be set to "BRI" or "BRIII", in order to
## automatically load them.
## treatLfc: treat log fold change for treat function when determining
## differentially expressed genes (with an FDR adjusted p.value cutoff of 0.01).
## adjMethod: p-value adjust method passed to p.adjust function to be applied to
## genes p-value calculation in order to define differentially expressed ones.
## pAdjCutOff: cutoff value to determine differentially expressed genes.
## ... : adittional parameters passed to dEnricher function.
IFA <- function(exprMatrix, classes, genesets=NULL, SEAcutoff=0.01, GSEAcutoff=0.01, br=NULL, treatLfc=0, adjMethod="fdr", pAdjCutOff=0.01, ...) {
    stopifnot(length(unique(classes)) == 2);
    
    class1 <- unique(classes)[[1]];
    class2 <- unique(classes)[[2]];
    
    # Delete genes which have less than 50% of lectures for each sample
    class1Data <- exprMatrix[, classes == class1];
    class2Data <- exprMatrix[, classes == class2];
    
    otherParams <- list(...);
    
    if (!missing(genesets)) {
	print("Using own gene sets");
	mgsz_gs <- genesets;
    } else if (all(c("ontology", "genome") %in% names(otherParams))) {
	print("Loading dEnricher gene sets");
	# Must use dEnricher data base
	genesets <- NULL;
	
	ontology <- otherParams$ontology;
	genome   <- otherParams$genome;
	
	if (!"RData.location" %in% names(otherParams)) {
	    # dEnricher default RData.location
	    RData.location <- as.list(args(dEnricher))$RData.location;
	} else {
	    RData.location <- otherParams$RData.location;
	}
	
	mgsz_gs <- dRDataLoader(paste("org.", genome, ".eg", ontology, sep=""), RData.location=RData.location, genome=genome, ontology=ontology)$gs;
	mgsz_gs <- lapply(mgsz_gs, unique);
    } else {
	print("Loading Gene Ontology gene sets");
	genesets <- loadGO();
	mgsz_gs <- genesets;
    }

    print("Starting SEA analysis");
    deRes <- DEnricher(exprMatrix, classes, genesets, treatLfc, br, SEAcutoff, adjMethod, pAdjCutOff, ...);
    colnames(deRes)[-1] <- paste("SEA", colnames(deRes)[-1], sep="_");
    
    print("Starting GSEA analysis");
    mgszRes <- Mgsz(exprMatrix, classes, mgsz_gs, GSEAcutoff);
    colnames(mgszRes)[-1] <- paste("GSEA", colnames(mgszRes)[-1], sep="_");
    
    ifaRes <- merge(deRes, mgszRes, by.x="setID", by.y="gene.sets", all=!F, suffixes=c("SEA","GSEA"));
    ifaRes$setID <- as.character(ifaRes$setID);
    ifaRes$Name <- unlist(lapply(ifaRes$setID, Term));
    return(ifaRes);
}

######## Usage Example:
# source("PaperCode/Pam50.R");
# 
# nki <- loadBCDataset(libname="nki", verbose=F);
# 
# # Use only genes with EntrezID
# nki <- nki[!is.na(nki$genes$EntrezGene.ID), ];
# 
# # Get subtype classification by PAM50
# nki <- PAM50Call(nki, std="median", verbose=F);
# 
# # EntrezID as rownames
# rownames(nki) <- nki$genes$EntrezGene.ID;
# nki <- nki[!is.na(rownames(nki)),];
# 
# basalData <- nki$M[, rownames(nki$subtype)[nki$subtype == "Basal"]];
# lumaData  <- nki$M[, rownames(nki$subtype)[nki$subtype == "LumA"]];
# exprMatrix <- cbind(basalData, lumaData);
# classes <- c(rep("Basal", ncol(basalData)), rep("LumA", ncol(lumaData)));
# table(classes);
# rm(nki); rm(basalData); rm(lumaData);
# 
# # this will test FA only for the first 75 gene sets from Gene Ontology
# ifaResults <- IFA(exprMatrix, classes, loadGO()[1:75], treatLfc=0.2);
