require(org.Hs.eg.db);

loadGO <- function() {
    go <- org.Hs.egGO2ALLEGS;
    go <- as.list(go);
    go <- lapply(go, unique);
    return(go);
}
