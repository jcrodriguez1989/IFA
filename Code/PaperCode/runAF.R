library(Biobase);
rm(list=ls());

#####
## Every timing benchmark was done over a single cores
## Intel(R) Xeon(R) CPU           X5650  @ 2.67GHz
#####

## Getting the data sets for analysis
source("Pam50.R");

getData <- function(libname, C1="Basal", C2="LumA") {
    myData <- loadBCDataset(libname=libname, verbose=F);
    
    # Use only genes with EntrezID
    myData <- myData[!is.na(myData$genes$EntrezGene.ID), ];
    
    # Get subtype classification by PAM50
    myData <- PAM50Call(myData, std="median", verbose=F);
    
    # EntrezID as rownames
    rownames(myData) <- myData$genes$EntrezGene.ID;
    myData <- myData[!is.na(rownames(myData)),];
    
    result <- list();
    C1data <- myData$M[, rownames(myData$subtype)[myData$subtype == C1]];
    C2data <- myData$M[, rownames(myData$subtype)[myData$subtype == C2]];
    
    
    # Saving genes which have at least 50% of lectures from the samples
    C1data <- C1data[apply(C1data, 1, function(x) sum(is.na(x))) < (ncol(C1data)/2),];
    C2data <- C2data[apply(C2data, 1, function(x) sum(is.na(x))) < (ncol(C2data)/2),];
    finalGenes <- intersect(rownames(C1data), rownames(C2data));
    
    print("****************************");
    print(libname);
    print(paste(ncol(C1data), C1, "subjects"));
    print(paste(ncol(C2data), C2, "subjects"));
    print(paste("BRIII", length(finalGenes)));
    
    result[[C1]] <- C1data[finalGenes,];
    result[[C2]] <- C2data[finalGenes,];
    
    return(result);
}

## It will print the number of subjects and BRIII of each data set (table 1)
libnames <-c ("vdx", "nki", "transbig", "upp", "unt", "mainz");
datasets <- lapply(libnames, getData);
names(datasets) <- libnames;

## Get differentially expressed genes for SEA analysis
getDiffGenes <- function(x, treatLfc=1) {
    C1 <- names(x)[[1]];
    C2 <- names(x)[[2]];
    
    stopifnot(rownames(x[[ C1 ]]) == rownames(x[[ C2 ]]));
    exprsX <- cbind(x[[ C1 ]], x[[ C2 ]]);
    
    # Gene to gene linear model using limma
    # Design matrix
    PAM50Call <- data.frame(PAM50.Call=c(rep(C1, ncol(x[[ C1 ]])), rep(C2, ncol(x[[ C2 ]]))));
    rownames(PAM50Call) <- c(colnames(x[[ C1 ]]), colnames(x[[ C2 ]]));
    design <- model.matrix(~PAM50.Call-1, data=PAM50Call);
    colnames(design) <- gsub("PAM50.Call", "", colnames(design));
    
    # Adjust the model
    fit <- lmFit(exprsX, design);
    
    # Specifying interest contrasts
    contrast.matrix <- makeContrasts(contrasts=paste(C1, C2, sep="-"), levels=design);
    
    # Ebayes correction
    fit <- treat(contrasts.fit(fit, contrast.matrix), lfc=treatLfc);
    # Adjusted pvalues
    fit$p.adjust <- apply(fit$p.value, 2, p.adjust, method="fdr");

    return(fit);
}

# With these treatLfc values, for each dataset I get around 5% of the genes
# (using a fixed fdr p.adjust <= 0.01)
# "vdx"    "nki"    "transbig"    "upp"    "unt"    "mainz" 
#  0.75     0.2         0.6        0.3     0.25      0.45
treatLfcs <- rbind(libnames, treatVal=c(.75, .2, .6, .3, .25, .45));

fits <- apply(treatLfcs, 2, function(x) getDiffGenes(datasets[[ x[1] ]], treatLfc=as.numeric(x[2]) ));
names(fits) <- libnames;

## It will print the number of DE genes for each data set and its percentage over BRIII (table 3)
invisible(lapply(names(fits), function(x) {
    actFit <- fits[[ x ]];
    
    total <- nrow(actFit);
    diff <- sum(actFit$p.adjust <= 0.01);
    percDiff <- diff * 100 / total;
    
    print(c(x, diff, total, round(percDiff,1)));
}))

deGenes <- lapply(fits, function(actualFit) { return(rownames(actualFit)[actualFit$p.adjust <= 0.01]) } );
names(deGenes) <- libnames;
deGenesInters <- apply(combn(libnames, 2), 2, function(x) {
    length(intersect(deGenes[[ x[[1]] ]], deGenes[[ x[[2]] ]]))
} );
names(deGenesInters) <- apply(combn(libnames, 2), 2, paste, collapse="_");

## It will print the intersection of DE genes between data sets (table 3)
deGenesInters;
rm(deGenesInters);

## total diff genes intersection (table 3)
length(Reduce(intersect, deGenes));

## total diff genes union (table 3)
length(Reduce(union, deGenes));
rm(deGenes);

## DAVID
library(RDAVIDWebService);

# briii: wether use bri or briii
runDAVID <- function(fit, email, briii=F) {
    pAdjCutOff <- 0.01;
    dif <- fit$p.adjust <= pAdjCutOff;
    dif <- unique(rownames(dif)[dif]);
    
    if (missing(email)) {
        stop("You should provide a registered DAVID email");
    }
    
    david <- DAVIDWebService$new(email=email, url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
    
    david$addList(dif, idType="ENTREZ_GENE_ID", listName="DEGenes", listType="Gene");
    
    setHttpProtocolVersion(david, version="HTTP/1.0");
    getHttpProtocolVersion(david);
    
    if (briii) {
	BRIII <- unique(rownames(fit)); length(BRIII)
	david$addList(BRIII, idType="ENTREZ_GENE_ID", listName="BRIII", listType="Background");
    }
    
    setTimeOut(david, 900000)
    getTimeOut(david)
    
    categories <- c("GOTERM_BP_ALL", "GOTERM_CC_ALL", "GOTERM_MF_ALL");
    david$setAnnotationCategories(categories);
    
    davidRes <- david$getFunctionalAnnotationChart(threshold=1, count=0L);
    
    enrFilterValue <- 0.05;
    davidRes$Category <- gsub("_.*", "", gsub("GOTERM_", "", davidRes$Category));
    davidRes$Enriched <- davidRes$PValue <= enrFilterValue;
    davidRes$GO_ID <- gsub("~.*", "", davidRes$Term);
    rownames(davidRes) <- davidRes$GO_ID;
    
    return(davidRes);
}

# this line takes 6.7 mins to execute.
WD_BRIII <- lapply(fits, function(x) runDAVID(x, briii=!F, email="PUT YOUR MAIL HERE!!"));
names(WD_BRIII) <- names(fits);
save(WD_BRIII, file="WD_BRIII.RData");
rm(WD_BRIII);

# this line takes 4.6 mins to execute.
WD_BRI <- lapply(fits, function(x) runDAVID(x, briii=F, email="PUT YOUR MAIL HERE!!"));
names(WD_BRI) <- names(fits);
save(WD_BRI, file="WD_BRI.RData");
rm(WD_BRI);


runRDavid <- function(fit, genesets, briii=F, testSide="greater") {
#     br <- unique(names(as.list(org.Hs.egGO)));
    br <- unique(unlist(genesets)); length(br)
    
    if (briii) {
	br <- intersect(br, unique(rownames(fit))); length(br)
    }
    
    pAdjCutOff <- 0.01;
    dif <- fit$p.adjust <= pAdjCutOff;
    dif <- unique(rownames(dif)[dif]); length(dif)
    dif <- intersect(dif, br); length(dif)
    
    pvals <- do.call(c, lapply(genesets, function(gset) {
	cont <- matrix(0, nrow=2, ncol=2);
	gset <- unique(intersect(gset, br));
	
	# see https://david.ncifcrf.gov/content.jsp?file=functional_annotation.html#fisher
	cont[1,1] <- length(intersect(dif, gset));
	cont[1,2] <- length(gset) - cont[1,1];
	cont[2,1] <- length(dif) - cont[1,1];
	cont[2,2] <- length(br) - sum(cont);
	
	cont[1,1] <- cont[1,1]-1 # EASE Score;
	
	if (cont[1,1] < 0) {
	    return(NA);
	} else {
	    return(fisher.test(cont, alternative=testSide)$p.value);
	}
    }));
    
    davidRes <- data.frame(GO_ID=names(genesets), pval=pvals);
    davidRes$p.adj <- p.adjust(pvals, "fdr");
    
    enrFilterValue <- 0.01;
    davidRes$Enriched <- davidRes$pval <= enrFilterValue;
    return(davidRes);
}

## In order to reproduce this publication, ensure you have org.Hs.eg.db v3.0.0
library(org.Hs.eg.db);
if (packageVersion("org.Hs.eg.db") != "3.0.0") {
    warning("In order to correctly reproduce this paper, get org.Hs.eg.db version 3.0.0");
}

loadGO <- function() {
    go <- org.Hs.egGO2ALLEGS
    go <- as.list(go)
    go <- lapply(go, unique);
    return(go);
}
go <- loadGO();

# this line takes 2.5 mins to execute.
RD_BRIII <- lapply(fits, function(x) runRDavid (x, go, briii=!F));
names(RD_BRIII) <- names(fits);
save(RD_BRIII, file="RD_BRIII.RData");
rm(RD_BRIII);

# this line takes 2.9 mins to execute.
RD_BRI <- lapply(fits, function(x) runRDavid (x, go, briii=F));
names(RD_BRI) <- names(fits);
save(RD_BRI, file="RD_BRI.RData");
rm(RD_BRI);

## Now we can execute "comparingCutoffs.R" in order to check the selected cutoff for RD


## GOstats
library(GOstats);

# https://github.com/Bioconductor-mirror/Category/blob/master/R/hyperGTest-methods.R
runGOstats <- function(fit, briii=F) {
    pAdjCutOff <- 0.01;
    dif <- fit$p.adjust <= pAdjCutOff;
    dif <- unique(rownames(fit)[dif]);

    params <- new("GOHyperGParams",
	geneIds=dif, # diff genes
	annotation="org.Hs.eg.db",
# 	pvalueCutoff=filterValue, # its the pval at which the term is considered enriched
	conditional=!FALSE
    );
    
    br <- "univ";
    if (briii) {
	params@universeGeneIds <- unique(rownames(fit));
    }

    ontologies <- c("BP", "CC", "MF");
    GOstats <- Reduce(rbind, lapply(ontologies, function(ont) {
	ontology(params) <- ont;

        testDirection(params) <- "over";
        resOver <- hyperGTest(params);
        filterValue <- resOver@pvalueCutoff;
        resOver <- summary(resOver, pvalue=1.1); # 1.1 in order to get all results
        colnames(resOver)[1] <- "GO_ID";
        resOver$Enriched <- resOver$Pvalue < filterValue;
        
	res <- resOver;
        res$ont <- ont;
        
	return(res);
    }));
    return(GOstats);
}

# this line takes 10.8 mins to execute.
GOstats_BRIII <- lapply(fits, function(x) runGOstats(x, brii=!F));
names(GOstats_BRIII) <- names(fits);
save(GOstats_BRIII, file="GOstats_BRIII.RData");
rm(GOstats_BRIII);

# this line takes 11.3 mins to execute.
GOstats_BRI <- lapply(fits, function(x) runGOstats(x, brii=F));
names(GOstats_BRI) <- names(fits);
save(GOstats_BRI, file="GOstats_BRI.RData");
rm(GOstats_BRI);


## dEnricher
library(dnet);

runDEnricher <- function(fit, method, test="HypergeoTest", verbose=F, filtering=F) {
    pAdjCutOff <- 0.01;
    dif <- fit$p.adjust <= pAdjCutOff;
    dif <- unique(rownames(fit)[dif]);
    
    sizeRange=c(10, 1000);
    min.overlap=3;
    if (!filtering) {
	sizeRange=c(0, 99999999);
	min.overlap=0;
    }
    
    # default parameters except for ontology algorithm used.
    # Note: HypergeoTest is default
    res_bp <- dEnricher(data=dif,
                identity="entrez",
                genome="Hs",
                ontology="GOBP",
                sizeRange=sizeRange,
                min.overlap=min.overlap,
                test=test,
                ontology.algorithm=method,
                verbose=verbose,
                p.adjust.method="BY");
    res_mf <- dEnricher(data=dif,
                identity="entrez",
                genome="Hs",
                ontology="GOMF",
                sizeRange=sizeRange,
                min.overlap=min.overlap,
                test=test,
                ontology.algorithm=method,
                verbose=verbose,
                p.adjust.method="BY");
    res_cc <- dEnricher(data=dif,
                identity="entrez",
                genome="Hs",
                ontology="GOCC",
                sizeRange=sizeRange,
                min.overlap=min.overlap,
                test=test,
                ontology.algorithm=method,
                verbose=verbose,
                p.adjust.method="BY");
    
    # enrichment cutoff is 0.01 as in dEnrichers Reference Manual it uses this
    # cutoff for elim.pvalue parameter
    enr_coff <- 0.01;
    res_bp <- data.frame(res_bp$set_info, pvalue=res_bp$pvalue, adjp=res_bp$adjp, Enriched=res_bp$pvalue < enr_coff);
    res_mf <- data.frame(res_mf$set_info, pvalue=res_mf$pvalue, adjp=res_mf$adjp, Enriched=res_mf$pvalue < enr_coff);
    res_cc <- data.frame(res_cc$set_info, pvalue=res_cc$pvalue, adjp=res_cc$adjp, Enriched=res_cc$pvalue < enr_coff);
    
    res <- rbind(res_bp, res_mf, res_cc);
    colnames(res)[[1]] <- "GO_ID";
    return(res);
}

## HyperGeo test
# dEnricher output
# 'org.Hs.eg' (from https://github.com/hfang-bristol/RDataCentre/blob/master/dnet/1.0.7/org.Hs.eg.RData?raw=true) has been loaded into the working environment (at 2016-06-09 12:37:19)
# 'org.Hs.egGO__' (from https://github.com/hfang-bristol/RDataCentre/blob/master/dnet/1.0.7/org.Hs.egGOBP.RData?raw=true) has been loaded into the working environment (at 2016-06-09 12:37:23)

# this line takes 5.817822 mins to execute.
dE_hyp_none <- lapply(fits, function(x) runDEnricher(x, method="none"));
names(dE_hyp_none) <- names(fits);
save(dE_hyp_none, file="dE_hyp_none.RData");
rm(dE_hyp_none);

# this line takes 30.7636 mins to execute.
dE_hyp_lea <- lapply(fits, function(x) runDEnricher(x, method="lea"));
names(dE_hyp_lea) <- names(fits);
save(dE_hyp_lea, file="dE_hyp_lea.RData");
rm(dE_hyp_lea);

# this line takes 33.56225 mins to execute.
dE_hyp_elim <- lapply(fits, function(x) runDEnricher(x, method="elim"));
names(dE_hyp_elim) <- names(fits);
save(dE_hyp_elim, file="dE_hyp_elim.RData");
rm(dE_hyp_elim);

# dE_hyp_pc <- lapply(fits, function(x) runDEnricher(x, method="pc"));
# names(dE_hyp_pc) <- names(fits);
# Could not test this algorithm because it has a bug with this output:
# Error in if (var.exp == 0) { : missing value where TRUE/FALSE needed


## binomial test
# this line takes 7.060218 mins to execute.
dE_bin_none <- lapply(fits, function(x) runDEnricher(x, method="none", test="BinomialTest"));
names(dE_bin_none) <- names(fits);
save(dE_bin_none, file="dE_bin_none.RData");
rm(dE_bin_none);

# this line takes 31.83038 mins to execute.
dE_bin_lea <- lapply(fits, function(x) runDEnricher(x, method="lea", test="BinomialTest"));
names(dE_bin_lea) <- names(fits);
save(dE_bin_lea, file="dE_bin_lea.RData");
rm(dE_bin_lea);

# this line takes 32.93243 mins to execute.
dE_bin_elim <- lapply(fits, function(x) runDEnricher(x, method="elim", test="BinomialTest"));
names(dE_bin_elim) <- names(fits);
save(dE_bin_elim, file="dE_bin_elim.RData");
rm(dE_bin_elim);

# dE_bin_pc <- lapply(fits, function(x) runDEnricher(x, method="pc", test="BinomialTest"));
# names(dE_bin_pc) <- names(fits);


## Fisher test
# this line takes 5.584001 mins to execute.
dE_fis_none <- lapply(fits, function(x) runDEnricher(x, method="none", test="FisherTest"));
names(dE_fis_none) <- names(fits);
save(dE_fis_none, file="dE_fis_none.RData");
rm(dE_fis_none);

# this line takes 28.90785 mins to execute.
dE_fis_lea <- lapply(fits, function(x) runDEnricher(x, method="lea", test="FisherTest"));
names(dE_fis_lea) <- names(fits);
save(dE_fis_lea, file="dE_fis_lea.RData");
rm(dE_fis_lea);

# this line takes 20.64207 mins to execute.
dE_fis_elim <- lapply(fits, function(x) runDEnricher(x, method="elim", test="FisherTest"));
names(dE_fis_elim) <- names(fits);
save(dE_fis_elim, file="dE_fis_elim.RData");
rm(dE_fis_elim);

# dE_fis_pc <- lapply(fits, function(x) runDEnricher(x, method="pc", test="FisherTest"));
# names(dE_fis_pc) <- names(fits);


## GSEA
calculatePvals <- function(dataset) {
    C1 <- names(dataset)[1];
    C2 <- names(dataset)[2];
    stopifnot(rownames(dataset[[C1]]) == rownames(dataset[[C2]]));
    exprsX <- cbind(dataset[[C1]], dataset[[C2]]);
    exprsX <- apply(exprsX, 2, as.numeric);
    
    # Gene to gene linear model
    # Design matrix
    PAM50Call <- data.frame(PAM50.Call=c(rep(C1, ncol(dataset[[C1]])), rep(C2, ncol(dataset[[C2]]))));
    rownames(PAM50Call) <- c(colnames(dataset[[C1]]), colnames(dataset[[C2]]));
    design <- model.matrix(~PAM50.Call-1, data=PAM50Call);
    colnames(design) <- gsub("PAM50.Call", "", colnames(design));
    
    # Adjust linear model
    fit <- lmFit(exprsX, design);
    
    # Specifying interest contrast
    contrast.matrix <- makeContrasts(contrasts=paste(C1, C2, sep="-"), levels=design);
    
    # Ebayes correction
    fit <- eBayes(contrasts.fit(fit, contrast.matrix));
    # adjust pvalues
    fit$p.adjust <- apply(fit$p.value, 2, p.adjust, method="fdr");

    return(fit);
}

writeCls <- function(dataset, savefile) {
    savefile <- paste(sep="", savefile, ".cls");
    conActual <- file(savefile, open="w");
    C1 <- names(dataset)[1];
    C2 <- names(dataset)[2];
    C1count <- ncol(dataset[[1]]);
    C2count <- ncol(dataset[[2]]);
    
    writeLines(paste(sep=" ", C1count + C2count, 2, 1), conActual);
    writeLines(paste(sep=" ", "#", C1, C2), conActual);
    writeLines(paste(paste(collapse=" ", rep(C1, C1count)), paste(collapse=" ", rep(C2, C2count))), conActual);
    close(conActual);
    
    return(savefile);
}

writeGct <- function(dataset, savefile) {
    savefile <- paste(sep="", savefile, ".gct");
    conActual <- file(savefile, open="w");

    stopifnot(nrow(dataset[[1]]) == nrow(dataset[[2]]));
    
    dataset[[1]][is.na(dataset[[1]])] <- as.character(NaN);
    dataset[[2]][is.na(dataset[[2]])] <- as.character(NaN);
    
    stopifnot(all(rownames(dataset[[1]]) == rownames(dataset[[2]])));
    writeLines(paste(sep="", "#1.2\n", nrow(dataset[[1]]), "\t", ncol(dataset[[1]]) + ncol(dataset[[2]])), conActual);
    close(conActual);
    
    bindedData <- cbind(dataset[[1]], dataset[[2]]);
    aa <- data.frame(rownames(bindedData), NA, bindedData);
    colnames(aa) <- c("NAME", "Description", colnames(bindedData));
    
    write.table(aa, append=T, sep="\t", quote=F, row.names=F, file=savefile);
    
    return(savefile);
}

writeGmt <- function(genesets, savefile) {
    savefile <- paste(sep="", savefile, ".gmt");
    conActual <- file(savefile, open="w");
    
    invisible(lapply(names(genesets), function(gset) {
	term <- as.vector(unlist(Term(gset)));
	genes <- genesets[[gset]];
	
	resString <- c(gset, term, genes);
	
	writeLines(resString, con=conActual, sep='\t');
	writeLines("", con=conActual);
    }));
    
    close(conActual);
    
    return(savefile);
}

writeRnk <- function(dataset, savefile, ranking="negLog") {
    fit <- calculatePvals(dataset);
    
    if (ranking=="negLog") { 
	fileName <- paste(savefile, "_logPval.rnk", sep="");
	toSave <- data.frame(-log(fit$p.adjust));
    } else if(ranking=="Pval") {
	fileName <- paste(savefile, "1_Pval.rnk", sep="");
	toSave <- data.frame(1-fit$p.adjust);
    } else if(ranking=="tscore") {
	fileName <- paste(savefile, "tscore.rnk", sep="");
	toSave <- data.frame(fit$t);
    } else {
	print("ranking must be one of negLog, Pval, tscore");
	return(0);
    }
    
    toSave <- toSave[order(toSave[,1], decreasing=T),, drop=F];
    
    write.table(toSave, file=fileName, quote=F, sep="\t", col.names=F);
    
    return(fileName);
}

runGsea <- function(gseaExec, dataset, genesets, preRank=F, genePerm=F, weight="1", ranking="negLog") {
    tmpDir <- tempfile(pattern=as.character(Sys.getpid()));
    dir.create(tmpDir);
    print(tmpDir);
    dir.create(paste(tmpDir, "resDir", sep=("/")));
    resDir <- paste(tmpDir, "resDir", sep=("/"));
    
    C1 <- names(dataset)[1];
    C2 <- names(dataset)[2];
    stopifnot(all(rownames(C1) == rownames(C2)));
    stopifnot(weight %in% c("0", "1", "2"));
    
    file.copy(gseaExec, paste(tmpDir, "/gsea.jar", sep=""));
    gseaExec <- paste(tmpDir, "/gsea.jar", sep="");
    
    exprMatrix <- cbind(dataset[[1]], dataset[[2]]);
    classes <- c(rep(C1, ncol(dataset[[1]])), rep(C2, ncol(dataset[[2]])));
    
    gmtFile <- writeGmt(genesets, paste(tmpDir, "genesets", sep="/"));
    
    weight <- switch(weight,
	"0"="classic",
	"1"="weighted",
	"2"="weighted_p2"
    )
    
    if (preRank) {
	outFile <- paste(tmpDir, paste(C1, C2, "preRank", ranking, weight, sep="_"), sep="/");
	
	rnkFile <- writeRnk(dataset, outFile, ranking); print("rnk Done");
	
	gseaCommand <- paste("java -cp ", gseaExec, " -Xmx4g xtools.gsea.GseaPreranked -rnk ", rnkFile, " -gmx ", gmtFile, " -collapse false -mode Max_probe -norm meandiv -nperm 1000 -scoring_scheme ", weight, " -rpt_label report -include_only_symbols true -make_sets true -plot_top_x 0 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ", resDir, " -gui false", sep="");
	
	C1 <- "pos";
	C2 <- "neg";
	
    } else {
	outFile <- paste(tmpDir, paste(C1, C2, ifelse(genePerm, "gP", "pP"), weight, sep="_"), sep="/");
	
	perm <- ifelse(genePerm, "gene_set", "phenotype");
	
	gctFile <- writeGct(dataset, outFile); print("gct Done");
	clsFile <- writeCls(dataset, outFile); print("cls Done");
	
	gseaCommand <- paste("java -cp ", gseaExec, " -Xmx4g xtools.gsea.Gsea -res ", gctFile, " -cls ", clsFile, "#", paste(C1, "_versus_", C2, sep=""), " -gmx ", gmtFile, " -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute ", perm, " -rnd_type no_balance -scoring_scheme ", weight, " -rpt_label report -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 0 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out ", resDir, " -gui false", sep="");
# 	 > /dev/null;
    }
    
    print("to execute");
    system(gseaCommand);
    print("done");
    C1File <- dir(resDir, pattern=paste("gsea", C1, ".xls", sep=".*"), recursive=!F, full.names=!F);
    C2File <- dir(resDir, pattern=paste("gsea", C2, ".xls", sep=".*"), recursive=!F, full.names=!F);

    C1res <- read.csv(file=C1File, sep="\t"); C1res <- cbind(C1res, REG=rep(C1, nrow(C1res)));
    C2res <- read.csv(file=C2File, sep="\t"); C2res <- cbind(C2res, REG=rep(C2, nrow(C2res)));

    gseaRes <- rbind(C1res, C2res);
    
    if (preRank | genePerm) {
	gseaRes$Enriched <- gseaRes$FDR.q.val < 0.05;
    } else {
	gseaRes$Enriched <- gseaRes$FDR.q.val < 0.25;
    }
    
    gmtUsed <- dir(resDir, pattern=".*gene_sets.gmt", recursive=!F, full.names=!F);
    attr(gseaRes, "gmtUsed") <- gmtUsed;
    
    return(gseaRes);
}

gseaExec <- "PUT GSEA .JAR EXECUTABLE FULL FILE PATH";
## In order to reproduce this paper, use GSEA v2-2.2.1
# for example gseaExec <- "~/gsea2-2.2.1.jar";

# this line takes 1.76 hours to execute.
SMgp0 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, genePerm=!F, weight="0"));
names(SMgp0) <- names(datasets);
save(SMgp0, file="SMgp0.RData");
# rm(SMgp0); # dont remove this because we are going to use it in mGSZ

# this line takes 1.84 hours to execute.
SMgp1 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, genePerm=!F, weight="1"));
names(SMgp1) <- names(datasets);
save(SMgp1, file="SMgp1.RData");
rm(SMgp1);

# this line takes 1.85 hours to execute.
SMgp2 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, genePerm=!F, weight="2"));
names(SMgp2) <- names(datasets);
save(SMgp2, file="SMgp2.RData");
rm(SMgp2);


# this line takes 1.72 hours to execute.
SMpp0 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, genePerm=F, weight="0"));
names(SMpp0) <- names(datasets);
save(SMpp0, file="SMpp0.RData");
rm(SMpp0);

# this line takes 1.84 hours to execute.
SMpp1 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, genePerm=F, weight="1"));
names(SMpp1) <- names(datasets);
save(SMpp1, file="SMpp1.RData");
rm(SMpp1);

# this line takes 1.84 hours to execute.
SMpp2 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, genePerm=F, weight="2"));
names(SMpp2) <- names(datasets);
save(SMpp2, file="SMpp2.RData");
rm(SMpp2);


# this line takes 1.52 hours to execute.
SMpr_log_p_0 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, preRank=!F, weight="0", ranking="negLog"));
names(SMpr_log_p_0) <- names(datasets);
save(SMpr_log_p_0, file="SMpr_log_p_0.RData");
rm(SMpr_log_p_0);

# this line takes 1.3 hours to execute.
SMpr_log_p_1 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, preRank=!F, weight="1", ranking="negLog"));
names(SMpr_log_p_1) <- names(datasets);
save(SMpr_log_p_1, file="SMpr_log_p_1.RData");
rm(SMpr_log_p_1);

# this line takes 1.3 hours to execute.
SMpr_log_p_2 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, preRank=!F, weight="2", ranking="negLog"));
names(SMpr_log_p_2) <- names(datasets);
save(SMpr_log_p_2, file="SMpr_log_p_2.RData");
rm(SMpr_log_p_2);


# this line takes 1.5 hours to execute.
SMpr_1_p_0 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, preRank=!F, weight="0", ranking="Pval"));
names(SMpr_1_p_0) <- names(datasets);
save(SMpr_1_p_0, file="SMpr_1_p_0.RData");
rm(SMpr_1_p_0);

# this line takes 1.34 hours to execute.
SMpr_1_p_1 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, preRank=!F, weight="1", ranking="Pval"));
names(SMpr_1_p_1) <- names(datasets);
save(SMpr_1_p_1, file="SMpr_1_p_1.RData");
rm(SMpr_1_p_1);

# this line takes 1.5 hours to execute.
SMpr_1_p_2 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, preRank=!F, weight="2", ranking="Pval"));
names(SMpr_1_p_2) <- names(datasets);
save(SMpr_1_p_2, file="SMpr_1_p_2.RData");
rm(SMpr_1_p_2);


# this line takes 1.5 hours to execute.
SMpr_tScore_0 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, preRank=!F, weight="0", ranking="tscore"));
names(SMpr_tScore_0) <- names(datasets);
save(SMpr_tScore_0, file="SMpr_tScore_0.RData");
rm(SMpr_tScore_0);

# this line takes 1.5 hours to execute.
SMpr_tScore_1 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, preRank=!F, weight="1", ranking="tscore"));
names(SMpr_tScore_1) <- names(datasets);
save(SMpr_tScore_1, file="SMpr_tScore_1.RData");
rm(SMpr_tScore_1);

# this line takes 1.5 hours to execute.
SMpr_tScore_2 <- lapply(datasets, function(dataset) runGsea(gseaExec, dataset, go, preRank=!F, weight="2", ranking="tscore"));
names(SMpr_tScore_2) <- names(datasets);
save(SMpr_tScore_2, file="SMpr_tScore_2.RData");
rm(SMpr_tScore_2);

## mGSZ
library(mGSZ);

runMGsz <- function(dataset, genesets) {
    C1 <- names(dataset)[1];
    C2 <- names(dataset)[2];
    stopifnot(all(rownames(C1) == rownames(C2)));
    
    exprMatrix <- cbind(dataset[[1]], dataset[[2]]);
    classes <- c(rep(C1, ncol(dataset[[1]])), rep(C2, ncol(dataset[[2]])));
    
    mGSZ.obj <- mGSZ(x=exprMatrix, y=genesets, l=classes);
    # in mGSZ min=19 => GS > 19, otherwise in GSEA min=19 => GS >= 19

    mGszRes <- mGSZ.obj$mGSZ;
    
    mGszRes$Enriched <- mGszRes$pvalue <= 0.01;
    
    return(mGszRes);
}

# this line takes 21.11 hours to execute.
mGSZ_5_inf <- lapply(datasets, function(x) runMGsz(x, go));
names(mGSZ_5_inf) <- names(fits);
save(mGSZ_5_inf, file="mGSZ_5_inf.RData");
rm(mGSZ_5_inf);

load.gmt.data <- function(gmt.file.path){
  tmp <- readLines(gmt.file.path)
  gsets <- list()
  for(i in 1:length(tmp)){
    t <- strsplit(tmp[i],'\t')[[1]]
    gsets[[t[1]]] <- t[3:length(t)]
  }
  return (lapply(gsets, unique))
}

# In order to check mGSZ with GSEAs 15-500 limits, we are going to use the genesets that GSEA used previously (after its filtering)
# this line takes 9.4 hours to execute.
mGSZ_15_500 <- lapply(names(datasets), function(x) {
    actGS <- load.gmt.data(attr(SMgp0[[x]], "gmtUsed"));
    runMGsz(datasets[[x]], actGS);
});
names(mGSZ_15_500) <- names(fits);
save(mGSZ_15_500, file="mGSZ_15_500.RData");
rm(mGSZ_15_500);

rm(go);


########### Lets run IFA for all datasets
setwd("../");
source("IFA.R");
setwd("PaperCode/");
load("tcgaInput.RData");

# these lines takes 22 hours to execute
gsets <- loadGO();
ifaRes <- lapply(names(datasets), function(dsetName) {
    print(paste("Starting:", dsetName));
    dset <- datasets[[dsetName]];
    exprMatrix <- cbind(dset[[1]], dset[[2]]);
    classes <- c(rep("Basal", ncol(dset[[1]])), rep("LumA", ncol(dset[[2]])));
    treatLfc <- as.numeric(treatLfcs[2, treatLfcs[1,] == dsetName ]);
    res <- IFA(exprMatrix, classes, gsets, treatLfc=treatLfc);
});
names(ifaRes) <- names(datasets);

ifaRes$tcga <- IFA(tcga$M, tcga$labels, gsets, treatLfc=1);
save(ifaRes, file="ifaRes.RData");
rm(ifaRes);

########### Lets run IFA for prostate cancer datasets
allGenes <- unique(unlist(gsets)); length(allGenes);
pCancerDsets <- list();

# Camcap
library("prostateCancerCamcap");
actData <- camcap;
exprMatrix <- exprs(camcap);

classes <- as.character(pData(phenoData(actData))$Sample_Group);
classes <- gsub("CRPC", "other", gsub("Tumour", "other", classes));

rownames(exprMatrix) <- as.character(fData(actData)$Entrez_Gene_ID);
exprMatrix <- exprMatrix[rownames(exprMatrix) %in% allGenes, ];
exprMatrix <- avereps(exprMatrix);
exprMatrix <- exprMatrix[rowSums(is.na(exprMatrix)) < ncol(exprMatrix)/2, ];

pCancerDsets$Camcap <- list(exprMatrix=exprMatrix, classes=classes, treatLfc=0.2);
rm(exprMatrix); rm(classes); rm(actData);

# taylor
library("prostateCancerTaylor");
actData <- taylor;
exprMatrix <- exprs(actData);

classes <- as.character(pData(phenoData(actData))$Sample_Group);

exprMatrix <- log2(exprMatrix);
matches <- match(rownames(exprMatrix), aux$ID);
exprMatrix <- exprMatrix[!is.na(matches),];
matches <- matches[!is.na(matches)];
rownames(exprMatrix) <- aux[matches, "Entrez"];

aux <- classes %in% c("normal adjacent benign prostate", "prostate cancer");
exprMatrix <- exprMatrix[,aux];
classes <- classes[aux];

exprMatrix <- exprMatrix[rownames(exprMatrix) %in% allGenes, ];
exprMatrix <- avereps(exprMatrix);
exprMatrix <- exprMatrix[rowSums(is.na(exprMatrix)) < ncol(exprMatrix)/2, ];
classes <- gsub("normal adjacent benign prostate", "benign", gsub("prostate cancer", "other", classes));

pCancerDsets$taylor <- list(exprMatrix=exprMatrix, classes=classes, treatLfc=0.15);
rm(aux); rm(matches); rm(exprMatrix); rm(classes); rm(actData);

# Varambally
library("prostateCancerVarambally");

actData <- varambally;
exprMatrix <- exprs(actData);

classes <- as.character(pData(phenoData(actData))$Sample_Group);
classes <- gsub("Metastatic", "other", gsub("Tumour", "other", classes));

exprMatrix <- log2(exprMatrix);
rownames(exprMatrix) <- as.character(fData(actData)$ENTREZ_GENE_ID);
exprMatrix <- exprMatrix[rownames(exprMatrix) %in% allGenes, ];
exprMatrix <- avereps(exprMatrix);
exprMatrix <- exprMatrix[rowSums(is.na(exprMatrix)) < ncol(exprMatrix)/2, ];

# get_treatVal(exprMatrix, classes, adjMethod="fdr"); # it doesnt get any enrichment
pCancerDsets$Varambally <- list(exprMatrix=exprMatrix, classes=classes, treatLfc=0.2, adjMethod="none");
rm(exprMatrix); rm(classes); rm(actData);

# Grasso
library("prostateCancerGrasso");
actData <- grasso;
exprMatrix <- exprs(actData);

classes <- as.character(pData(phenoData(actData))$Group);
classes <- gsub("CastrateResistant", "other", gsub("HormomeDependant", "other", classes));

rownames(exprMatrix) <- as.character(fData(actData)$GENE);
exprMatrix <- exprMatrix[rownames(exprMatrix) %in% allGenes, ];
exprMatrix <- avereps(exprMatrix);
exprMatrix <- exprMatrix[rowSums(is.na(exprMatrix)) < ncol(exprMatrix)/2, ];

pCancerDsets$Grasso <- list(exprMatrix=exprMatrix, classes=classes, treatLfc=0.45);
rm(exprMatrix); rm(classes); rm(actData); rm(allGenes);

# these lines takes 15.88674 hours to execute
pCancerIfaRes <- lapply(names(pCancerDsets), function(dsetName) {
    print(paste("Starting:", dsetName));
    dset <- pCancerDsets[[dsetName]];
    exprMatrix <- dset$exprMatrix;
    classes <- dset$classes;
    treatLfc <- dset$treatLfc;
    
    if ("adjMethod" %in% names(dset)) {
	res <- IFA(exprMatrix, classes, gsets, treatLfc=treatLfc, adjMethod=dset$adjMethod);
    } else {
	res <- IFA(exprMatrix, classes, gsets, treatLfc=treatLfc);
    }
    return(res);
});
names(pCancerIfaRes) <- names(pCancerDsets);
save(pCancerIfaRes, file="pCancerIfaRes.RData");

rm(gsets);

########### check for correct saving
rm(list=ls());

resFiles <- c(
"WD_BRI.RData", 
"WD_BRIII.RData", 
"RD_BRI.RData", 
"RD_BRIII.RData", 
"GOstats_BRI.RData", 
"GOstats_BRIII.RData", 
"dE_hyp_none.RData", 
"dE_hyp_lea.RData", 
"dE_hyp_elim.RData", 
"dE_bin_none.RData", 
"dE_bin_lea.RData", 
"dE_bin_elim.RData", 
"dE_fis_none.RData", 
"dE_fis_lea.RData", 
"dE_fis_elim.RData", 
"SMgp0.RData", 
"SMgp1.RData", 
"SMgp2.RData", 
"SMpp0.RData", 
"SMpp1.RData", 
"SMpp2.RData", 
"SMpr_log_p_0.RData", 
"SMpr_log_p_1.RData", 
"SMpr_log_p_2.RData", 
"SMpr_1_p_0.RData", 
"SMpr_1_p_1.RData", 
"SMpr_1_p_2.RData", 
"SMpr_tScore_0.RData", 
"SMpr_tScore_1.RData", 
"SMpr_tScore_2.RData", 
"mGSZ_5_inf.RData", 
"mGSZ_15_500.RData");

mergeResults <- function(resFiles) {
    # checking some stuff
    invisible(lapply(resFiles, function(resFile) {
	print(resFile);
	loadedData <- get(load(resFile));
	stopifnot(length(loadedData) == 6); #check all datasets have results
	print(names(loadedData));
	print(lapply(loadedData, dim));
	lapply(loadedData, function(x) stopifnot(all(c("GO_ID", "Enriched") %in% colnames(x)) | all(c("NAME", "Enriched") %in% colnames(x)) | all(c("gene.sets", "Enriched") %in% colnames(x))));
    }))

    aux <- lapply(resFiles, function(resFile) {
	loadedData <- get(load(resFile));
	method <- gsub(".RData", "", resFile);
	
	go_id <- intersect(c("GO_ID", "NAME", "gene.sets"), colnames(loadedData[[1]]));
	loadedData <- lapply(loadedData, function(x) x[, c(go_id, "Enriched")]);
	resDframe <- suppressWarnings(Reduce(function(...) merge(..., by=go_id, all=!F), loadedData));
	colnames(resDframe) <- c("GO_ID", paste(method, names(loadedData), sep=":"));
	
	return(resDframe);
    });
    resDframe <- suppressWarnings(Reduce(function(...) merge(..., by="GO_ID", all=!F), aux));
    rm(aux);
    return(resDframe);
}

resDframe <- mergeResults(resFiles);
save(resDframe, file="allResults.RData");
## Now we can continue with file resAnalysis.R


