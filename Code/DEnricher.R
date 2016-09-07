require(dnet);

getDiffGenes <- function(exprMatrix, classes, treatLfc, adjMethod) {
    class1 <- unique(classes)[[1]];
    class2 <- unique(classes)[[2]];
    
    # Gene to gene linear model using limma
    # Design matrix
    dMatrix <- data.frame(d.Matrix=classes);
    rownames(dMatrix) <- colnames(exprMatrix);
    design <- model.matrix(~d.Matrix-1, data=dMatrix);
    colnames(design) <- gsub("d.Matrix", "", colnames(design));
    
    # Adjust the model
    fit <- lmFit(exprMatrix, design);
    
    # Specifying interest contrasts
    contrast.matrix <- makeContrasts(contrasts=paste(class1, class2, sep="-"), levels=design);
    
    # Ebayes correction
    fit <- treat(contrasts.fit(fit, contrast.matrix), lfc=treatLfc);
    # Adjusted pvalues
    fit$p.adjust <- apply(fit$p.value, 2, p.adjust, method=adjMethod);

    return(fit);
}

get_treatVal <- function(exprMatrix, classes, dePerc=5, presicion=0.05, adjMethod="fdr", pAdjCutOff=0.01) {
    if (dePerc > 1) dePerc <- dePerc/100;
    treats <- seq(0, 1, by=presicion);
    deCounts <- do.call(c, lapply(treats, function(pres) {
	fit <- getDiffGenes(exprMatrix, classes, pres, adjMethod);
	dif <- fit$p.adjust <= pAdjCutOff;
	actDePerc <- sum(dif) / length(dif);
	print(paste("Act treat:", pres, "- DE percentage:", round(actDePerc,6)));
	return(actDePerc);
    }));
    return(treats[[which.min(abs(deCounts - dePerc))]]);
}


runDEnricher <- function(dif, genesets, br, test, SEAcutoff) {
    GeneID <- dif
    
    genes.group <- GeneID[!is.na(GeneID)]
    set_filtered <- names(genesets);

    gs <- genesets;
    nSet <- length(gs)
    
    doFisherTest <- function(genes.group, genes.term, genes.universe) {
        genes.hit <- intersect(genes.group, genes.term)
        X <- length(genes.hit)
        K <- length(genes.group)
        M <- length(genes.term)
        N <- length(genes.universe)
        cTab <- matrix(c(X, K - X, M - X, N - M - K + X), nrow = 2, 
            dimnames = list(c("anno", "notAnno"), c("group", 
                "notGroup")))
        p.value <- ifelse(all(cTab == 0), 1, stats::fisher.test(cTab, 
            alternative = "greater")$p.value)
        return(p.value)
    }
    doHypergeoTest <- function(genes.group, genes.term, genes.universe) {
        genes.hit <- intersect(genes.group, genes.term)
        X <- length(genes.hit)
        K <- length(genes.group)
        M <- length(genes.term)
        N <- length(genes.universe)
        x <- X
        m <- M
        n <- N - M
        k <- K
        p.value <- ifelse(m == 0 || k == 0, 1, stats::phyper(x, 
            m, n, k, lower.tail = F, log.p = F))
        return(p.value)
    }
    doBinomialTest <- function(genes.group, genes.term, genes.universe) {
        genes.hit <- intersect(genes.group, genes.term)
        X <- length(genes.hit)
        K <- length(genes.group)
        M <- length(genes.term)
        N <- length(genes.universe)
        p.value <- ifelse(K == 0 || M == 0 || N == 0, 1, stats::pbinom(X, 
            K, M/N, lower.tail = F, log.p = F))
        return(p.value)
    }

    terms <- names(gs);
#     genes.universe <- as.numeric(unique(unlist(gs[terms])))
    genes.universe <- br;
    genes.group <- intersect(genes.universe, genes.group)
    if (length(genes.group) == 0) {
        warnings("There is no gene being used.\n")
        return(F)
    }
    
    pvals <- sapply(terms, function(term) {
	genes.term <- as.numeric(unique(unlist(gs[term])))
	p.value <- switch(test, FisherTest = doFisherTest(genes.group, 
	    genes.term, genes.universe), HypergeoTest = doHypergeoTest(genes.group, 
	    genes.term, genes.universe), BinomialTest = doBinomialTest(genes.group, 
	    genes.term, genes.universe))
    })
    
    pvals <- sapply(pvals, function(x) min(x, 1));
    pvals <- signif(pvals, digits = 2);
    
    seaRes <- data.frame(
	    setID=as.character(names(pvals)),
	    pval=as.numeric(pvals),
	    Enriched=pvals <= SEAcutoff);
    return(seaRes);
}

DEnricher <- function(exprMatrix, classes, genesets, treatLfc, br, SEAcutoff, adjMethod, pAdjCutOff, ...) {
    fit <- getDiffGenes(exprMatrix, classes, treatLfc, adjMethod);
    dif <- fit$p.adjust <= pAdjCutOff;
    dif <- unique(rownames(dif)[dif]); length(dif)
    print(paste("DE genes", length(dif), "of a total of", nrow(exprMatrix), "(", round(length(dif)/nrow(exprMatrix)*100,2), "%)"));

    otherParams <- list(...);
    
    if (!"test" %in% names(otherParams)) test <- "FisherTest";
    if (!"sizeRange" %in% names(otherParams)) sizeRange <- c(0, 999999);
    if (!"min.overlap" %in% names(otherParams)) min.overlap <- 0;

    if (is.null(genesets)) {
	pvals <- dEnricher(dif, test=test, sizeRange=sizeRange, min.overlap=min.overlap, ...);

	if (!"set_info" %in% names(pvals)) stop("Error in dEnricher.");
	
	seaRes <- data.frame(
	    setID=as.character(pvals$set_info[, "setID"]),
	    pval=as.numeric(pvals$pvalue),
	    p.adj=p.adjust(pvals$pvalue, "fdr"),
	    Enriched=pvals$pvalue <= SEAcutoff);
    } else {
        if (is.null(br)) {
	    br <- unique(unlist(genesets));
	    print(paste("Using BRI:", length(br), "genes."));
        } else if (length(br) == 1 & br == "I") {
	    br <- unique(unlist(genesets));
	    print(paste("Using BRI:", length(br), "genes."));
        } else if (length(br) == 1 & br == "III") {
	    br <- unique(rownames(exprMatrix));
	    print(paste("Using BRIII:", length(br), "genes."));
	} else if (length(br) > 1) {
	    br <- intersect(br, unique(unlist(genesets)));
	    print(paste("Using user provided BR:", length(br), "genes."));
	} else {
	    stop("Incorrect br option.");
	}
	
	dif <- intersect(br, dif);
	genesets <- lapply(genesets, intersect, br);
	genesets <- genesets[ !is.null(genesets) ];
	
	seaRes <- runDEnricher(dif, genesets, br, test=test, SEAcutoff=SEAcutoff);
    }
    
    print(paste(sum(seaRes$Enriched, na.rm=!F), " enriched terms", sep=""));
    return(seaRes);
}
