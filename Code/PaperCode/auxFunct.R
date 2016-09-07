# checks how contained are methods into other (inter-methods)
getContention <- function(resDframe, methName1, methName2) {
    dsets <- c("vdx","nki","transbig","upp","unt","mainz");
    contains <- do.call(c,lapply(dsets, function(dset) {
	meth1 <- resDframe[,paste(methName1, ":", dset, sep="")];
	meth2 <- resDframe[,paste(methName2, ":", dset, sep="")];
	meth1[is.na(meth1)] <- F;
	meth2[is.na(meth2)] <- F;
	
	print("-------------------------------");
	print(dset);
	res <- table(meth1, meth2);
	print(res);
	
        # method 2 enriched terms contained in method 1
        return(res[2,2] / sum(res[,2]));
    }));
    return(contains);
}

# Converts results list. Each list member is a database represented as a data.frame with genes x methods (logical value wether enriched or not).
resultsAsDataFrame <- function(results) {
    stopifnot(length(results) > 1);
    
    aa <- lapply(names(results), function(x) {
	colnames(results[[x]])[-1] <- paste(x, colnames(results[[x]])[-1], sep="");
	return(results[[x]])
    } );
    names(aa) <- names(results);
    aa <- Reduce(function(...) merge(..., by="NAME", all=!F), aa);
    
    return(aa);
}

# Given a method name and a data.frame, per gene it sums the number of databases in which it was enriched. If all values are NA then it returns NA
getSums <- function(method, dframe) {
    regexp <- paste(method, ":", sep=""); # paste(method, "$", sep="");
    actData <- dframe[,grep(regexp, colnames(dframe))];
    stopifnot(ncol(actData) == 6);
    
    result <- rowSums(actData, na.rm=!F);
    result[apply(actData, 1, function(y) all(is.na(y)))] <- NA;
    return(result);
}

library(GO.db);
# Internal Function
getDepth <- function (term, minDepth, allParents) {
    fstTerms <- c("GO:0008150", "GO:0003674", "GO:0005575"); # bp, mf, cc
    actualDepth <- 1;
    actualParents <- unlist(allParents[term]);
    
    if (term %in% fstTerms) {
	return(0);
    } else if (is.null(actualParents)) {
	return(NA);
    }
    
    if (minDepth) {
	while ((!fstTerms %in% actualParents) && (actualParents != "all")) {
	    actualDepth <- actualDepth +1;
	    actualParents <- unlist(allParents[unlist(actualParents)]);
	}
    } else {
	while (!all(actualParents %in% as.list(c("all", fstTerms)))) {
	    actualDepth <- actualDepth +1;
	    actualParents <- unlist(allParents[unlist(actualParents)]);
	}
    }
    
    return(actualDepth);
}

# Given a list of terms (GO IDs), it returns the depths of them.
# If minDepth is true then the shortest path to the root is calculated, otherwise, the longest path
getDepths <- function (terms, minDepth=TRUE) {
    allParents <- as.list(GOBPPARENTS);
    allParents <- c(allParents, as.list(GOMFPARENTS));
    allParents <- c(allParents, as.list(GOCCPARENTS));
    
    return(lapply(terms, function(x) { getDepth(x, minDepth, allParents) } ));
}

# gomatrix: matrix NxP (N: terms count, P: databases count) with values [0, 1]
# flevel: percentage [0,1] of enriched terms beneath the databases to include in heatmap
# penrich: percentage [0,1] of enriched terms permitted by database
# ... extra parameters for heatmap.2
library(gplots);
library(vegan);
filter.map <- function(gomatrix, flevel=0, penrich=0, methods=c(), dendrogram, col.dist="jaccard", row.dist=col.dist, ...) {
    if (missing(dendrogram)) dendrogram <- "column";
    # percentage per row of enrichment
    
    col.colors <- rep("white", ncol(gomatrix));
    if (length(methods) > 0) {
	colors <- rainbow(length(methods));
	colorIndexes <- do.call(rbind, lapply(1:length(methods), function(x) {
	    actMethod <- methods[[x]];
# 	    regexp <- paste(actMethod, "$", sep="");
	    regexp <- actMethod;
	    colIndexes <- grep(regexp, colnames(gomatrix));
	    cbind(colIndexes, rep(colors[[x]], length(colIndexes)));
	}));
	
	col.colors[ as.numeric(colorIndexes[,1]) ] <- colorIndexes[,2];
    }

    totalP <- rowMeans(gomatrix, na.rm=!F);
    # number of enriched terms per database (column)
    totalE <- colSums(gomatrix, na.rm=!F) / nrow(gomatrix);
    # delete terms under the desired cutoff
    go.idx <- which(totalP >= flevel);
    print(length(go.idx));
    enr.idx <- which(totalE >= penrich);
    print(length(enr.idx));
    
    # jaccard clustering per column (database)
    dd.s <- suppressWarnings(vegdist(t(gomatrix[go.idx, enr.idx]), col.dist, na.rm=!F));
    h.s <- hclust(dd.s, method="average");
    # jaccard clustering per row (term)
    ddr.s <- suppressWarnings(vegdist(gomatrix[go.idx,enr.idx], row.dist, na.rm=!F));
    ddr.s[is.na(ddr.s)] <- 0;
    hr.s <- hclust(ddr.s,method="average")

    gomatrix[is.na(gomatrix)] <- -1;
    colorsMap <- c("white", "orange1", "red1");
    
    heatmap.2(as.matrix(gomatrix[go.idx,enr.idx]), Rowv=as.dendrogram(hr.s), dendrogram=dendrogram, ColSideColors=col.colors[enr.idx], trace="none", labRow=rep("", nrow(gomatrix)), Colv=as.dendrogram(h.s), colsep=1:(ncol(gomatrix)-1), sepwidth=c(0.025, 0.025), key=F, col=colorsMap, breaks=length(colorsMap)+1, ...);
    
    # return selected terms and the clustering
    return(list(row.index = go.idx,row.dendrogram = as.dendrogram(hr.s),col.index=enr.idx,col.dendrogram=as.dendrogram(h.s)))
}


# returns te number of terms which return the queried disease
queryCtdbase <- function(terms, disease="breast") {
    if (length(terms)==0) return(FALSE);
    if (is.na(terms)) return(FALSE);

    urlfst <- "http://ctdbase.org/tools/batchQuery.go?inputType=go&inputTerms=";
    urllst <- "|mercury&report=diseases_inferred&format=tsv";
    
    print("Starting queries");
    
    inDisease <- do.call(c, lapply(terms, function(term) {
	aa <- tryCatch({read.csv(paste(urlfst, term, urllst, sep=""), sep="\t")[,-7]},
	error={function(x) return(NA)},
	warning={function(x) return(NA)});
	
	if (is.na(aa)) {
	    print(paste("NA returned for ", term, sep=""));
	    return(FALSE);
	}
	wDisease <- any(aa[,3] == term & grepl(disease, aa[, 4], ignore.case=!F));
	
	print(paste("***** Disease in term ", term, ": ", wDisease, sep=""));
	
	if (wDisease) {
	    print(unique(as.character(aa[grep(disease, aa[, 4], ignore.case=!F), 4])));
	}
	
	return(wDisease);
    }));
    return(inDisease);
}

# GO Semantic Similarity
library(GOSemSim);

# measures <- c("Resnik", "Lin", "Rel", "Jiang", "Wang");
getGOSim <- function(measure) {
    print(paste("Get", measure, "measures."));
    onts <- c("BP", "CC", "MF");
    GOSim <- lapply(onts, function(actOnt) {
	print(paste("Ontology:", actOnt));
	actResDframe <- resDframe[Ontology(as.character(resDframe[,1])) == actOnt,];
	combs <- combn(colnames(actResDframe)[-1], 2);
	aa <- apply(combs, 2, function(actComb) {
	    comb1aux <- actResDframe[, actComb[[1]] ];
	    comb1aux[ is.na(comb1aux) ] <- FALSE;
	    comb1GO <- as.character(actResDframe[ comb1aux, 1 ]);
	    
	    comb2aux <- actResDframe[, actComb[[2]] ];
	    comb2aux[ is.na(comb2aux) ] <- FALSE;
	    comb2GO <- as.character(actResDframe[ comb2aux, 1 ]);
	    
	    res <- mgoSim(comb1GO, comb2GO, ont=actOnt, organism="human", measure=measure);
	    return(res);
	});
	dists <- matrix(0, nrow=ncol(actResDframe)-1, ncol=ncol(actResDframe)-1);
	colnames(dists) <- colnames(actResDframe)[-1];
	rownames(dists) <- colnames(actResDframe)[-1];
	
	for (i in 1:length(aa)) {
	    actComb <- combs[,i];
	    dists[ actComb[[1]], actComb[[2]] ] <- aa[[ i ]];
	    dists[ actComb[[2]], actComb[[1]] ] <- aa[[ i ]];
	}
	dists <- as.dist(dists);
	return(dists);
    })
    names(GOSim) <- onts;
    return(GOSim);
}

mergeIfaRes <- function(ifaResList) {
    aux <- lapply(names(ifaResList), function(actDbase) {
	actRes <- ifaResList[[actDbase]];
	actRes <- actRes[,c("setID", "GSEA_Enriched", "SEA_Enriched")];
	colnames(actRes) <- c("setID", paste(actDbase, c("GSEA", "SEA"), sep=":"));
	
	return(actRes);
    });
    resDframe <- suppressWarnings(Reduce(function(...) merge(..., by="setID", all=!F), aux));
    rm(aux);
    return(resDframe);
}

getWhoEnrich <- function(ifaResDFrame) {
    dsets <- unique(gsub(":.*", "", colnames(ifaResDFrame)[-1]));
    whoEnrich <- do.call(cbind, lapply(dsets, function(actDset) {
	gseaData <- ifaResDFrame[,paste(actDset, "GSEA", sep=":")] == 1;
	seaData  <- ifaResDFrame[,paste(actDset, "SEA", sep=":")] == 1;
	
	res <- do.call(c,lapply(1:length(gseaData), function(i) {
	    actRes <- c(gseaData[[i]], seaData[[i]]);
	    if (sum(is.na(actRes)) == 2) {
		return("NA");
	    } else if (sum(is.na(actRes)) == 1) {
		actRes[is.na(actRes)] <- F;
	    }
	    if (all(actRes)) {
		return("both");
	    } else if (actRes[[1]]) {
		return("GSEA");
	    } else if (actRes[[2]]) {
		return("SEA");
	    } else {
		return("None");
	    }
	}))
	return(res);
    }))
    colnames(whoEnrich) <- dsets;
    rownames(whoEnrich) <- ifaResDFrame[,1];
    
    return(whoEnrich);
}

ifaHeatmap <- function(enrichDFrame) {
    # jaccard clustering per column (database)
    dd.s <- suppressWarnings(vegdist(t(enrichDFrame), "jaccard", na.rm=!F));
    h.s <- hclust(dd.s, method="average");
    # jaccard clustering per row (term)
    ddr.s <- suppressWarnings(vegdist(enrichDFrame, "jaccard", na.rm=!F));
    ddr.s[is.na(ddr.s)] <- 0;
    hr.s <- hclust(ddr.s,method="average");
    
    colorsMap <- c("orange1", "red1");
    if (any(is.na(enrichDFrame))) {
	warning("White color represents NA.");
	colorsMap <- c("white", "orange1", "red1");
    }
    
    heatmap.2(as.matrix(enrichDFrame), Rowv=as.dendrogram(hr.s), dendrogram="column", trace="none", labRow=rep("", nrow(enrichDFrame)),Colv=as.dendrogram(h.s), colsep=1:(ncol(enrichDFrame)-1), sepwidth=c(0.025, 0.025), key=F, col=colorsMap, breaks=length(colorsMap)+1);
    return(NA);
}

