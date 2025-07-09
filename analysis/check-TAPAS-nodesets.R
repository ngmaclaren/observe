library(optNS)

dl.dw <- readRDS("../data/datalist-TAPAS-fmri-doublewell.rds")
dl.sis <- readRDS("../data/datalist-TAPAS-fmri-SIS.rds")
dl.gr <- readRDS("../data/datalist-TAPAS-fmri-genereg.rds")
dl.wc <- readRDS("../data/datalist-TAPAS-fmri-wilsoncowan.rds")
tsdir <- "/projects/academic/naokimas/neil/brains-ns50/"
tsfiles <- list.files(path = tsdir, pattern = ".txt")
pIDs <- gsub(".txt", "", tsfiles)

whichpID <- 2
pID <- pIDs[whichpID]

df <- read.table(paste0(tsdir, tsfiles[whichpID]), sep = " ")
nr <- nrow(df)
X.train <- as.matrix(df[seq_len(nr-100), ])
X.test <- as.matrix(df[seq(nr-100+1, nr, by = 1), ])

ns.dw <- dl.dw[[whichpID]]
ns.sis <- dl.sis[[whichpID]]
ns.gr <- dl.gr[[whichpID]]
ns.wc <- dl.wc[[whichpID]]

                                        # plot the time series
palette("Paired")
pdf("../img/tstest.pdf")
par(mar = c(5, 5, 1, 1), pty = "s")
                                        # This is the time series. It's constant for all conditions
simtime <- as.numeric(rownames(X.test))
matplot(
    type = "l",
    x = simtime, xlab = "t",
    y = X.test, ylab = "x",
    main = pID,
    cex.lab = 1.5, cex.axis = 1.25,
    lty = 1, lwd = 0.25, col = adjustcolor("black", 0.1)
)
lines(simtime, rowMeans(X.test), lty = 1, lwd = 6, col = "black")
                                        # doublewell
matlines(
    x = simtime,
    y = X.test[, ns.dw$opts[[1]]$vs],
    lty = 1, lwd = 3, col = 1
)
lines(simtime, rowMeans(X.test[, ns.dw$opts[[1]]$vs]), lty = 1, lwd = 6, col = 2)
                                        # genereg
matlines(
    x = simtime,
    y = X.test[, ns.gr$opts[[1]]$vs],
    lty = 1, lwd = 3, col = 3
)
lines(simtime, rowMeans(X.test[, ns.gr$opts[[1]]$vs]), lty = 1, lwd = 6, col = 4)
                                        # legend
legend(
    "topleft", bty = "n", cex = 1.25, lwd = 3, col = c(2, 4),
    legend = c("Double-well", "Gene-regulatory")
)
dev.off()

                                        # what about the relative errors on the orignial dynamics?
sapply(ns.dw[c("opts", "rands")], function(x) summary(get_error(x)))
sapply(ns.sis[c("opts", "rands")], function(x) summary(get_error(x)))
sapply(ns.gr[c("opts", "rands")], function(x) summary(get_error(x)))
sapply(ns.wc[c("opts", "rands")], function(x) summary(get_error(x)))
