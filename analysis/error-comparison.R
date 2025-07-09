library(optNS)
library(sfsmisc)

save_plots <- TRUE # FALSE
weighted <- FALSE # TRUE
directed <- TRUE # FALSE

dynamics <- "doublewell" # mutualistic SIS genereg

networks <- c( # unweighted, undirected
    "dolphin", "proximity", "celegans", "euroroad", "email", 
    "drosophila", "reactome", "route_views", "spanish", "foldoc", "tree_of_life", "word_assoc",
    "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr" 
)
labels <- c(
    "Dolphin", "Proximity", "Metabolic", "Road", "Email", "FlyBi", "Reactome", "Route views",
    "Spanish words", "FOLDOC", "Tree of life", "English words", "Enron", "Marker Cafe", "Prosper",
    "ER", "BA", "HK", "GKK", "LFR"
)
imgfile <- paste0("../img/error-comparison/", dynamics, ".pdf")

if(weighted) {
    networks <- c(
        "windsurfers", "macaques", "train_terrorists", "highschool", "drug", "residence_hall", "netsci_weighted",
        "proximity_weighted", "gap_junction_herm", "intl_trade"
    )
    labels <- c(
        "Windsurfer contact", "Macaque dominance", "Terrorist contact", "High school friendship",
        "Drug interaction", "University friendship", "Coauthorship", "Weighted proximity", "Neuronal", "Export"
    )
    imgfile <- paste0("../img/error-comparison/", dynamics, "-weighted.pdf")
}
if(directed) {
    networks <- c(
        "canton", "physician_trust", "email_company", "flamingo", "ecoli", "yeast", "usair", "jung-c", "email_uni",
        "faa"
    )
    labels <- c(
        "Trophic", "Physicians", "Email (manufacturer)", "Flamingo", "Escherichia coli", "Saccharomyces cerevisiae",
        "Flights", "JUNG", "Email (university)", "Air route"
    )
    imgfile <- paste0("../img/error-comparison/", dynamics, "-directed.pdf")
}

allns <- lapply(
    networks,
    function(network) readRDS(paste0("../data/ns-", network, "_", dynamics, ".rds"))
)
allerrors <- lapply(allns, function(ns) lapply(ns, get_error))

makedf <- function(dl, network) { # use mapply
    df <- data.frame(opt = dl$opt, fixed = dl$fixed, rand = dl$rand)
    rdf <- reshape(
        df,
        varying = c("opt", "fixed", "rand"),
        v.names = "error",
        timevar = "ns.type",
        times = c("opt", "fixed", "rand"),
        new.row.names = seq_len(dim(df)[1]*dim(df)[2]), direction = "long"
    )
    rdf <- rdf[, -which(colnames(rdf) == "id")]
    rdf$network <- network
    rdf$ns.type <- factor(rdf$ns.type, levels = c("rand", "opt", "fixed")) # order supports plotting
    return(rdf)
}

errordfs <- mapply(makedf, allerrors, networks, SIMPLIFY = FALSE)
## df <- do.call(rbind, errordfs) ## don't think I want to do this

## df <- errordfs[[1]] # for testing
## dev.new(height = 1.25, wd = 7)

plotit <- function(df, label, showmean = FALSE) {
    plot(
        df$error, jitter(rep(1, nrow(df)), amount = 0.05),
        col = as.numeric(df$ns.type), pch = as.numeric(df$ns.type) - 1, log = "x",
        axes = FALSE, xlab = "error", ylab = "", ylim = 1 + c(-0.075, 0.075)
    )
    if(showmean) points(mean(df$error[df$ns.type == "rand"]), 1, pch = 3, cex = 2)
    eaxis(1, n.axp = 1, cex.axis = 0.95)
    axis(2, at = 1, labels = label, tick = FALSE, las = 2)
}

ht <- 8
wd <- 7
mfrow <- c(20, 1)

if(weighted | directed) {
    ht <- 4
    mfrow <- c(10, 1)
}

palette(c("#33d17a", "#3584e4", "#ff7800"))
if(save_plots) {
    ## pdf(paste0("../img/error-comparison/", dynamics, "-wmean.pdf"), height = ht, width = wd)
    pdf(imgfile, height = ht, width = wd)
} else {
    dev.new(width = wd, height = ht)
}
par(mfrow = mfrow, mar = c(1.5, 7, 0, 0.5), mgp = c(3, 0.5, 0))
mapply(plotit, errordfs, labels)
if(save_plots) dev.off()
