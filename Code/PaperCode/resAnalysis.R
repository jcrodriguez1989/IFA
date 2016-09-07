rm(list=ls());
library(Cairo);
library(ggplot2);
library(reshape2);

source("auxFunct.R");
load("allResults.RData");

# enriched terms by dataset/method
methods <- unique(gsub(":.*", "", colnames(resDframe[,-1])));
enrichedbyDset <- do.call(cbind, lapply(methods, function(actMethod) {
    cols <- grep(paste("^", actMethod, ":.*", sep=""), colnames(resDframe));
    
    stopifnot(length(cols) == 6);
    aa <- colSums(resDframe[,cols], na.rm=!F);
    names(aa) <- gsub(".*:", "", names(aa));
    stopifnot(names(aa) == c("vdx","nki","transbig","upp","unt","mainz"));
    return(aa);
}))
colnames(enrichedbyDset) <- methods;

methodsMedian <- median(enrichedbyDset[enrichedbyDset > 0]); methodsMedian

# this is the order I want for boxplot
boxPlotOrder <- c("SMgp0","SMgp1","SMgp2","SMpp0","SMpp1","SMpp2","SMpr_tScore_0","SMpr_tScore_1","SMpr_tScore_2","SMpr_log_p_0","SMpr_log_p_1","SMpr_log_p_2","SMpr_1_p_0","SMpr_1_p_1","SMpr_1_p_2","mGSZ_15_500","mGSZ_5_inf","RD_BRI","RD_BRIII","WD_BRI","WD_BRIII","GOstats_BRI","GOstats_BRIII","dE_hyp_none","dE_hyp_lea","dE_hyp_elim","dE_bin_none","dE_bin_lea","dE_bin_elim","dE_fis_none","dE_fis_lea","dE_fis_elim");

stopifnot(length(boxPlotOrder) ==  ncol(enrichedbyDset));
enrichedbyDset <- enrichedbyDset[,boxPlotOrder];
rm(boxPlotOrder);

aux <- apply(resDframe[,-1],2, function(x) sum(!is.na(x)));
aux2 <- do.call(rbind,lapply(colnames(enrichedbyDset), function(actMethod) {
    aa <- do.call(rbind,lapply(rownames(enrichedbyDset), function(actDS) {
	actComb <- paste(actMethod, actDS, sep=":");
	c(actMethod, actDS, enrichedbyDset[actDS, actMethod], aux[[actComb]]);
    }))
    actQuant <- quantile(as.numeric(aa[, 3]), probs=c(.25,.5,.75));
    bb <- c(actQuant, IQR=actQuant[3]-actQuant[1], Var=var(as.numeric(aa[, 3])), Mean=mean(as.numeric(aa[, 3])));
    print(actMethod)
    print(round(bb, 2));
    return(aa);
}))
colnames(aux2) <- c("Method", "Data Set", "Enriched", "Analyzed");

## It will print for each method the number of enriched and analyzed terms (table S3)
aux2;
# write.table(aux2, file="../../Table_S3.csv", sep="\t", row.names=F, quote=F);
rm(aux); rm(aux2);

enrichedbyDset <- melt(enrichedbyDset);
colnames(enrichedbyDset) <- c("dataset", "method", "count");
# giving color to each method
colors <- rainbow(9);
fills <- c(rep(colors[[1]], 3), rep(colors[[2]], 12), rep(colors[[3]], 2), rep(colors[[4]], 2), rep(colors[[5]], 2), rep(colors[[6]], 2), rep(colors[[7]],3), rep(colors[[8]],3), rep(colors[[9]],3));

p <- ggplot();
p <- p + geom_boxplot(data=enrichedbyDset, aes(x=method, y=count), fill=fills);
p <- p + theme_bw();
p <- p + theme(axis.text.x=element_text(angle=45, hjust=1));
p <- p + xlab("Method") + ylab("Number of enriched terms");
p <- p + theme(panel.grid.major=element_blank()); # white background
p <- p + theme(panel.grid.minor=element_blank()); # white background
p <- p + geom_hline(aes(yintercept=methodsMedian));
rm(colors); rm(fills);

# check that columns are correct
# p <- p + scale_x_discrete(labels=c("SMgp0", "SMgp1", "SMgp2", "SMpp0", "SMpp1", "SMpp2", "SMpr tScore 0", "SMpr tScore 1", "SMpr tScore 2", "SMpr -log(p) 0", "SMpr -log(p) 1", "SMpr -log(p) 2", "SMpr 1-p 0", "SMpr 1-p 1", "SMpr 1-p 2", "mGSZ [15,500)", "mGSZ [5,∞)", "RD BRI", "RD BRIII", "WD BRI", "WD BRIII", "GOstats BRI", "GOstats BRIII", "dE_H_none", "dE_H_lea", "dE_H_elim", "dE_b_none", "dE_b_lea", "dE_b_elim", "dE_F_none", "dE_F_lea", "dE_F_elim"));
# ggsave(p, file="../../../FIG2.PDF", device=cairo_pdf);
# ggsave(p, file="../../../FIG2.JPEG");
rm(methodsMedian);
rm(p);

# check that results are contained into other (inter-methods): SEAs and mGSZ
mean(getContention(resDframe, "WD_BRI", "WD_BRIII")); # 0.9909482
mean(getContention(resDframe, "RD_BRI", "RD_BRIII")); # 0.9981257
mean(getContention(resDframe, "GOstats_BRI", "GOstats_BRIII")); # 0.8587952
mean(getContention(resDframe, "mGSZ_5_inf", "mGSZ_15_500")); # 0.9254076

# lets analyze dEnricher results
dEnrRes <- resDframe[, c(1, grep("dE_", colnames(resDframe)))];

mean(getContention(dEnrRes, "dE_fis_none", "dE_fis_lea"));  # 1
mean(getContention(dEnrRes, "dE_fis_none", "dE_fis_elim")); # 0.9929938

methods <- unique(gsub(":.*", "", colnames(dEnrRes)[-1]));
enrMethod <- do.call(cbind, lapply(methods, function(x) {
    colSums(dEnrRes[, grep(x, colnames(dEnrRes))], na.rm=!F);
}));

methods <- methods[grep("fis", methods)];
aux <- do.call(rbind, apply(combn(methods, 2), 2, function(actComb) {
    meth1 <- actComb[[1]];
    meth2 <- actComb[[2]];
    
    dsets <- c("vdx", "nki", "transbig", "upp", "unt", "mainz");
    meth1GS <- unique(do.call(c, lapply(dsets, function(x) {
	as.character(dEnrRes[ dEnrRes[, paste(meth1, x, sep=":")] & !dEnrRes[, paste(meth2, x, sep=":")], "GO_ID"]);
    })));
    
    meth2GS <- unique(do.call(c, lapply(dsets, function(x) {
	as.character(dEnrRes[ !dEnrRes[, paste(meth1, x, sep=":")] & dEnrRes[, paste(meth2, x, sep=":")], "GO_ID"]);
    })));
    
    meth1GS <- meth1GS[!is.na(meth1GS)];
    meth2GS <- meth2GS[!is.na(meth2GS)];
    
    res <- data.frame(
	total=c(length(meth1GS),
		length(meth2GS)),
	inCancer=c( sum(queryCtdbase(meth1GS), na.rm=!F),
		    sum(queryCtdbase(meth2GS), na.rm=!F))
    );
    rownames(res) <-  c(paste("in", meth1, "not", meth2, sep="_"), 
			paste("in", meth2, "not", meth1, sep="_"))
    
    return(res);
})); aux
#                    total inCancer
# in none not lea     39       34
# in lea not none      0        0
# in none not elim    52       44
# in elim not none     6        5
# in lea not elim     25       21
# in elim not lea     10        7
rm(methods); rm(enrMethod); rm(dEnrRes); rm(aux);

# densities with unified datasets
# methods filtering
resDframe <- resDframe[,-grep("SMgp0:|SMgp2:|SMpp0:|SMpp1:|SMpr_|mGSZ_15_500:|dE_hyp|dE_bin|dE_fis_lea|dE_fis_elim", colnames(resDframe))];
resDframe[,-1] <- apply(resDframe[,-1], 2, as.numeric);

resultsSumbyMethod <- resDframe[,1,drop=F];
impMethods <- unique(gsub(":.*", "", colnames(resDframe)))[-1]; impMethods;
resultsSumbyMethod <- cbind(resultsSumbyMethod, do.call(cbind, lapply(impMethods, function(method) getSums(method, resDframe))));
colnames(resultsSumbyMethod)[-1] <- impMethods;

# terms filtering
resultsSumbyMethod <- resultsSumbyMethod[rowSums(resultsSumbyMethod[,-1], na.rm=!F) > 0,];
resultsSumbyMethod$depths <- unlist(getDepths(as.character(resultsSumbyMethod$GO_ID)));
densities <- data.frame(c(by(resultsSumbyMethod[,-c(1,ncol(resultsSumbyMethod))], resultsSumbyMethod$depths, colSums, na.rm=!F)));
## It will print each methods number of enriched terms by depth (table S4)
colnames(densities) <- sort(unique(resultsSumbyMethod$depths)); densities

# in percentage
# densitiesPerc <- t(apply(densities, 1, function(x) (x*100)/sum(x)));
densitiesPerc <- densities;
densitiesPerc <- do.call(cbind, lapply(list(1:3, 4:5, 6:7, 8:ncol(densitiesPerc)), function(x) {
    rowSums(densitiesPerc[,x]);
}));
# densitiesPerc <- densitiesPerc / 6;

densitiesPerc <- 100*(densitiesPerc / rowSums(densitiesPerc)); densitiesPerc
# rowSums(densities); rowSums(densities)/6;

meltedDensities <- melt(data.frame(densitiesPerc, rownames(densitiesPerc)));
colnames(meltedDensities) <- c("Method", "Depth", "Percentage");

meltedDensities$Method <- gsub("mGSZ_5_inf", "mGSZ [5,∞)", gsub("_BRI", " BRI", gsub("dE_fis_none", "dE_F_none", as.character(meltedDensities$Method))));

p <- ggplot();
p <- p + geom_bar(data=meltedDensities, aes(x=Method, fill=as.factor(Depth), weight=Percentage), color="black") + coord_flip() + labs(list(y="Enrichment percentage", fill="GO tree\ndepth"));
p <- p + theme(panel.background=element_rect(fill="white")); # white background
# check depth intervals
# p <- p + scale_fill_manual(values=c("white", "peachpuff", "lightsalmon2", "indianred4"), labels=c("0-2", "3-4", "5-6", "7-10"));
# ggsave(p, file="../../../FIG3.PDF", device=cairo_pdf);
# ggsave(p, file="../../../FIG3.JPEG");
rm(p);

# Heatmap.
# methods filtering
impMethods <- c("SMpp2", "SMgp1", "mGSZ_5_inf", "WD_BRIII", "GOstats_BRIII", "RD_BRIII", "dE_fis_none");
regexp <- paste(impMethods, collapse="|");
resDframe <- cbind(GO_ID=resDframe$GO_ID, resDframe[,grep(regexp, colnames(resDframe))]);
# terms filtering
resDframe <- resDframe[rowSums(resDframe[,-1], na.rm=!F) > 0,];

# fixing some colnames for heatmap
resDframe2 <- resDframe;
impMethods2 <- unique(gsub(":.*", "", colnames(resDframe2)))[-1];

# jpeg("../../../FIG4_orig.JPG", height=1024, width=1024, quality=100);
# tiff("../../../FIG4_orig.TIFF", height=1024, width=1024, res=200);
# postscript("../../../FIG4_orig.EPS", height=1024, width=1024);
invisible(filter.map(resDframe2[,-1], 0, 0, impMethods2));
# dev.off();
rm(resDframe2); rm(impMethods2);

resultsSumbyMethod <- resDframe[,1,drop=F];
resultsSumbyMethod <- cbind(resultsSumbyMethod, do.call(cbind, lapply(impMethods, getSums, resDframe)));
colnames(resultsSumbyMethod)[-1] <- impMethods;
resultsSumbyMethod[,-1] <- resultsSumbyMethod[,-1]*100/6;

# 1:6*100 / 6
# [1]  16.66667  33.33333  50.00000  66.66667  83.33333 100.00000
percIn  <- 80;
percOut <- 20;

# Positive enrichment stability
sort(apply(resultsSumbyMethod[,-1], 2, function(x) round(sum(x > percIn, na.rm=!F) / sum(x > percOut, na.rm=!F),2)));

# Negative enrichment stability
sort(apply(resultsSumbyMethod[,-1], 2, function(x) round(sum(x < percOut, na.rm=!F) / sum(x < percIn, na.rm=!F),2)));

resultsSNa <- resultsSumbyMethod;
resultsSNa[is.na(resultsSNa)] <- 0;

methods <- colnames(resultsSNa)[-1];

# concordant temrs enriched by each method (i.e. in > 80 of dsets).
concterms <- lapply(methods, function(actualMethod) {
    aux <- resultsSNa[ resultsSNa[,actualMethod] > percIn, "GO_ID" ];
    return(as.character(aux));
});
names(concterms) <- methods;

# lets see that mGSZ tends to enrich almost all the EC over each other
mgsz <- "mGSZ_5_inf";
aa <- do.call(c,lapply(setdiff(methods, mgsz), function(x) {
    return(length(intersect(concterms[[x]], concterms[[mgsz]])) / length(concterms[[x]]));
}));
cbind(aa, setdiff(methods, mgsz))
mean(aa);
rm(mgsz); rm(aa);

DE <- "dE_fis_none";
aa <- do.call(c,lapply(setdiff(methods, DE), function(x) {
    return(length(intersect(concterms[[x]], concterms[[DE]])) / length(concterms[[x]]));
}));
cbind(aa, setdiff(methods, DE))
mean(aa);
rm(DE); rm(aa);

nonconcterms <- lapply(methods, function(actualMethod) {
    aux <- resultsSNa[ resultsSNa[,actualMethod] < percOut, "GO_ID" ];
    return(as.character(aux));
});
names(nonconcterms) <- methods;

getExclterms <- function(method, otherMethods) {
    aa <- Reduce(intersect, lapply(otherMethods, function(oMeth) {
	return(nonconcterms[[oMeth]]);
    }));
    return(intersect(concterms[[method]], aa));
}

# we delete WD for analysis hereafter, as RD replaces it
methods <- setdiff(methods, "WD_BRIII");
exclterms <- lapply(methods, function(actMethod) {
    otherMethods <- setdiff(methods, actMethod);
    print(actMethod);
    exclterms <- getExclterms(actMethod, otherMethods);
    print(length(exclterms));
    return(exclterms);
});
names(exclterms) <- methods;

# exclusively enriched terms
termsdepths <- data.frame(GO_ID=as.character(resultsSumbyMethod$GO_ID), depths=unlist(getDepths(as.character(resultsSumbyMethod$GO_ID))));

lapply(exclterms, function(x) { table(termsdepths[ termsdepths$GO_ID %in% x, "depths" ]) });

## this line requires internet connection in order to query the Comparative Toxigenomics Database
cancerInfo <- lapply(exclterms, queryCtdbase);

excltable <- rbind(do.call(c, lapply(exclterms, length)), do.call(c, lapply(cancerInfo, sum)));
excltable <- rbind(excltable, excltable[2,] / excltable[1,]); excltable

allInfo <- lapply(methods, function(actualMethod) {
    aux <- data.frame(cbind(exclterms[[actualMethod]], cancerInfo[[actualMethod]]));
    colnames(aux) <- c("GO_ID", "inCancer");
    aux <- merge(aux, termsdepths, by="GO_ID", all.x=!F);
})
names(allInfo) <- methods;

## It will print for each method its EET and if they are breast cancer related (table 4)
summarytable <- do.call(rbind, lapply(names(allInfo), function(actMethod) {
    actualMethod <- allInfo[[actMethod]];
    intervals <- list(0:2, 3:4, 5:6, 7:10, 11:50);
    actInfo <- do.call(cbind, lapply(intervals, function(actInt) {
	if (!any(actInt %in% actualMethod$depths)) {
	    return(c(0, 0));
	}
	aux <- actualMethod[ actualMethod$depths %in% actInt, ];
	return(c(nrow(aux), sum(aux$inCancer != FALSE)));
    }));
    rownames(actInfo) <- paste(actMethod, c("eet", "inCancer"), sep="_");
    colnames(actInfo) <- do.call(c, lapply(intervals, function(x) paste(min(x), max(x), sep="-")));
    return(actInfo);
}));
summarytable <- cbind(summarytable, total=rowSums(summarytable)); summarytable
#                        0-2 3-4 5-6 7-10 11-50 total
# SMpp2_eet                0   0   2    1     0     3
# SMpp2_inCancer           0   0   0    1     0     1
# SMgp1_eet                2  15   7    0     0    24
# SMgp1_inCancer           1   9   4    0     0    14
# mGSZ_5_inf_eet           0  26  27    8     0    61
# mGSZ_5_inf_inCancer      0  13  12    3     0    28
# GOstats_BRIII_eet        0   1   0    0     0     1
# GOstats_BRIII_inCancer   0   0   0    0     0     0
# RD_BRIII_eet             4   2   0    0     0     6
# RD_BRIII_inCancer        0   1   0    0     0     1
# dE_fis_none_eet          1   4   5    3     0    13
# dE_fis_none_inCancer     1   3   3    1     0     8


library(GO.db);
goterms <- Term(GOTERM);

finaltable <- do.call(rbind, lapply(names(allInfo), function(x) {
    aa <- allInfo[[x]];
    aa$Method <- x;
    aa$Ont <- do.call(c, lapply(as.character(aa$GO_ID), function(y) Ontology(y)));
    aa$Name <- do.call(c, lapply(as.character(aa$GO_ID), function(y) goterms[[y]]));
    aa <- aa[order(aa$depths),];
    return(aa);
}));
rm(goterms);

finaltable <- finaltable[, c("Method", "GO_ID", "inCancer", "depths", "Ont", "Name")];
# write.table(finaltable, file="../../Table_S5.csv", quote=F, row.names=F, sep="\t");

## IFA results analysis
load("ifaRes.RData");

ifaRes <- mergeIfaRes(ifaRes);
dim(ifaRes);
# terms filtering
ifaRes <- ifaRes[rowSums(ifaRes[,-1], na.rm=!F) > 0,];

whoEnrich <- getWhoEnrich(ifaRes);

aux <- whoEnrich;
aux[aux == "NA"] <- NA;
aux <- aux != "None";
aux <- data.frame(setID=as.character(ifaRes$setID), aux);
aux$depths <- unlist(getDepths(as.character(aux$setID)));

goterms <- Term(GOTERM);
aux$Ont <- do.call(c, lapply(as.character(aux$setID), function(y) Ontology(y)));
aux$Name <- do.call(c, lapply(as.character(aux$setID), function(y) goterms[y]));
rm(goterms);

anyEnrich <- whoEnrich;
anyEnrich[anyEnrich == "NA"] <- NA;
anyEnrich[anyEnrich == "None"] <- 0;
anyEnrich[anyEnrich != "0"] <- 1;
anyEnrich <- apply(anyEnrich,2,as.numeric);
rownames(anyEnrich) <- ifaRes$setID;

testData <- anyEnrich[,"tcga",drop=F];
trainingData <- anyEnrich[,!grepl("tcga", colnames(anyEnrich))];
# just keep gene sets enriched by tcga
trainingData <- trainingData[rowSums(testData, na.rm=!F) == 1,,drop=F];
nrow(trainingData);

sum(rowSums(trainingData, na.rm=!F) == 6) / nrow(trainingData);
sum(rowSums(trainingData, na.rm=!F) == 6);
# [1] 0.3325123
# [1] 270

sum(rowSums(trainingData, na.rm=!F) > 4) / nrow(trainingData); # EC
sum(rowSums(trainingData, na.rm=!F) > 4);
# [1] 0.4334975
# [1] 352

sum(rowSums(trainingData, na.rm=!F) == 0) / nrow(trainingData);
sum(rowSums(trainingData, na.rm=!F) == 0);
# [1] 0.1514778
# [1] 123

aux$label <- "";
aux[ aux$setID %in% names(which(rowSums(trainingData, na.rm=!F) > 4)), "label"] <- "CONCORDANT";
aux[ aux$setID %in% names(which(rowSums(trainingData, na.rm=!F) == 0)), "label"] <- "TCGA-EXCLUSIVE";
# write.table(aux, file="../../Table_S6.csv", quote=F, row.names=F, sep="\t");
rm(aux);

# mean(do.call(c,lapply(1:(ncol(anyEnrich)-1), function(x) {
#     act   <- anyEnrich[,x,drop=F];
#     other <- anyEnrich[,-x];
#     excl <- sum(rowSums(act, na.rm=!F) == 1 & rowSums(other, na.rm=!F) == 0);
#     excl / sum(rowSums(act, na.rm=!F) == 1);
# })))
# # [1] 0.1012566

# check if inCancer
concterms <- rownames(anyEnrich)[rowSums(anyEnrich) >= .8*ncol(anyEnrich)]; length(concterms);
inCancer <- queryCtdbase(concterms);
length(concterms); sum(inCancer, na.rm=!F); sum(inCancer, na.rm=!F) / length(concterms);
# [1] 445
# [1] 232
# [1] 0.5213483

exclterms <- do.call(rbind, lapply(colnames(anyEnrich), function(actDS) {
    actDSenr <- rownames(anyEnrich)[anyEnrich[, actDS] == 1];
    aux <- rownames(anyEnrich)[rowSums(anyEnrich, na.rm=!F) == 1];
    actDSexcl <- intersect(actDSenr, aux);
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    print(actDS); print(length(actDSexcl));
    inCancer <- queryCtdbase(actDSexcl);
    return(c(length(actDSexcl), sum(inCancer, na.rm=!F)));
}));
rownames(exclterms) <- colnames(anyEnrich);
exclterms <- cbind(exclterms, exclterms[,2] / exclterms[,1]); round(exclterms, 2)
#          [,1] [,2] [,3]
# vdx        75   38 0.51
# nki       170   84 0.49
# transbig  110   55 0.50
# upp        94   33 0.35
# unt        92   29 0.32
# mainz      58   30 0.52
# tcga      123   55 0.45

# jpeg("../../../FIG5_orig.JPG", height=1024, width=1024);
# tiff("../../../FIG5_orig.TIFF", height=1024, width=1024, res=200);
# postscript("../../../FIG5_orig.EPS", height=1024, width=1024);
ifaHeatmap(anyEnrich);
# dev.off();

## prostate Cancer IFA results
load("pCancerIfaRes.RData");

pCancerIfaRes <- mergeIfaRes(pCancerIfaRes);
dim(pCancerIfaRes);

# terms filtering
pCancerIfaRes <- pCancerIfaRes[rowSums(pCancerIfaRes[,-1], na.rm=!F) > 0,];

whoEnrich <- getWhoEnrich(pCancerIfaRes);

aux <- whoEnrich;
aux[aux == "NA"] <- NA;
aux <- aux != "None";
aux <- data.frame(setID=as.character(pCancerIfaRes$setID), aux);
aux$depths <- unlist(getDepths(as.character(aux$setID)));

goterms <- Term(GOTERM);
aux$Ont <- do.call(c, lapply(as.character(aux$setID), function(y) Ontology(y)));
aux$Name <- do.call(c, lapply(as.character(aux$setID), function(y) goterms[y]));
rm(goterms);
# write.table(aux, file="../../Table_S8.csv", quote=F, row.names=F, sep="\t");
rm(aux);

anyEnrich <- whoEnrich;
anyEnrich[anyEnrich == "NA"] <- NA;
anyEnrich[anyEnrich == "None"] <- 0;
anyEnrich[anyEnrich != "0"] <- 1;
anyEnrich <- apply(anyEnrich,2,as.numeric);
rownames(anyEnrich) <- pCancerIfaRes$setID;
rownames(whoEnrich) <- pCancerIfaRes$setID;

# check if inCancer
concterms <- rownames(anyEnrich)[rowSums(anyEnrich) >= .8*ncol(anyEnrich)];
inCancer <- queryCtdbase(concterms, disease="prostat");
length(concterms); sum(inCancer, na.rm=!F); sum(inCancer, na.rm=!F) / length(concterms);
# [1] 163
# [1] 99
# [1] 0.607362

exclterms <- do.call(rbind, lapply(colnames(anyEnrich), function(actDS) {
    actDSenr <- rownames(anyEnrich)[anyEnrich[, actDS] == 1];
    aux <- rownames(anyEnrich)[rowSums(anyEnrich, na.rm=!F) == 1];
    actDSexcl <- intersect(actDSenr, aux);
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
    print(actDS); print(length(actDSexcl));
    inCancer <- queryCtdbase(actDSexcl, disease="prostat");
    return(c(length(actDSexcl), sum(inCancer, na.rm=!F)));
}));
rownames(exclterms) <- colnames(anyEnrich);
exclterms <- cbind(exclterms, exclterms[,2] / exclterms[,1]); round(exclterms, 2)
#            [,1] [,2] [,3]
# Camcap      357  153 0.43
# taylor      448  194 0.43
# Varambally  212   98 0.46
# Grasso      263  112 0.43

# jpeg("../../../FIG6_orig.JPG", height=1024, width=1024);
# tiff("../../../FIG6_orig.TIFF", height=1024, width=1024, res=200);
# postscript("../../../FIG6_orig.EPS", height=1024, width=1024);
ifaHeatmap(anyEnrich);
# dev.off();

# pCancerterms <- rownames(anyEnrich)[rowSums(anyEnrich, na.rm=!F) == ncol(anyEnrich)];
pCancerterms <- rownames(anyEnrich)[rowSums(anyEnrich, na.rm=!F) > 2];
ctdbaseRes <- queryCtdbase(pCancerterms, disease="prostat");
sum(ctdbaseRes, na.rm=!F) / length(pCancerterms);
# [1] 0.6196319 in 4 datasets
# [1] 0.5978022 in at least 3 datasets
