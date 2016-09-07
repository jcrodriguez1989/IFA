library(vegan);

load("WD_BRI.RData");
load("RD_BRI.RData");
names(WD_BRI) <- paste("bri", names(WD_BRI), sep="_");
names(RD_BRI) <- paste("bri", names(RD_BRI), sep="_");
WD <- WD_BRI;
RD <- RD_BRI;
rm(WD_BRI); rm(RD_BRI);

load("WD_BRIII.RData");
load("RD_BRIII.RData");
names(WD_BRIII) <- paste("briii", names(WD_BRIII), sep="_");
names(RD_BRIII) <- paste("briii", names(RD_BRIII), sep="_");
WD <- c(WD, WD_BRIII);
RD <- c(RD, RD_BRIII);
rm(WD_BRIII); rm(RD_BRIII);

stopifnot(all(names(WD) == names(RD)));
pAdjusted <- F; # bool determining wether use padj or not for RD

## It will print the number of analyzed terms by WD, RD and their intersection (table S1)
allDists <- do.call(rbind, lapply(names(WD), function(dset) {
    actWD <- WD[[dset]];
    
    actRD <- RD[[dset]];
    actRD <- actRD[ !is.na(actRD$pval), ];
    
    commonGO <- intersect(actRD$GO_ID, actWD$GO_ID);
    
    print("*******************************");
    print(dset);
    print(paste("terms analyzed by RD:", nrow(actRD)));
    print(paste("terms analyzed by WD:", nrow(actWD)));
    print(paste("terms analyzed by both:", length(commonGO)));
    
    actRD <- actRD[ commonGO, ]; # filter the gsets analyzed by both methods
    actWD <- actWD[ commonGO, ];
    stopifnot(all(rownames(actRD) == rownames(actWD)));
    
    cOffs <- seq(0.005, 0.1, by=0.005);
    dists <- do.call(c, lapply(cOffs, function(cOff) {
	if (pAdjusted) {
	    actRD$Enriched <- actRD$p.adj <= cOff;
	} else {
	    actRD$Enriched <- actRD$pval <= cOff;
	}
# 	print(sum(actRD$Enriched==actWD$Enriched));
	enrichments <- rbind(as.numeric(actRD$Enriched), as.numeric(actWD$Enriched));
 	actDist <- vegdist(enrichments, "jaccard", na.rm=!F);
 	return(actDist);
    }));
    names(dists) <- cOffs;
    return(dists);
}));
rownames(allDists) <- names(WD);

## It will print Jaccard distances between RD and WD for every data set and cutoff value (table S2)
round(t(allDists),2);
round(colMeans(allDists),2);
which.min(colMeans(allDists));
# 0.01

rm(allDists); rm(pAdjusted);
rm(RD); rm(WD);
