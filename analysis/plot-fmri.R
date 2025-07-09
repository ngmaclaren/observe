dynamics <- "genereg"
TAPAS <- TRUE

if(TAPAS) {
    imgfile <- paste0("../img/epsilon-TAPAS-dist-", dynamics, ".pdf")
    dl <- readRDS(paste0("../data/datalist-TAPAS-fmri-", dynamics, ".rds"))
    flag <- readRDS("./TAPAS-flag4.rds")
} else {
    imgfile <- paste0("../img/epsilon-dist-", dynamics, ".pdf")
    dl <- readRDS(paste0("../data/datalist-fmri-", dynamics, ".rds"))
    flag <- sapply(dl, `[[`, "flag")
}
                                        # check this for each dynamics, then remove.
                                        # should be the same?
## length(which(sapply(dl, `[[`, "flag") == "use")) # 900

dlr <- dl[which(flag == "use")]

get_testerrors <- function(datalist) {
    errors <- datalist$testerror
    sapply(errors, mean)
}

error <- data.frame(t(sapply(dlr, get_testerrors)))

## Histogram version
optcolor <- adjustcolor("#3585e4")  
fixedcolor <- adjustcolor("#ff7800")
randcolor <- adjustcolor("#33d17a") 
frac <- TRUE # FALSE

hists <- list()
hists$all <- hist(unlist(error), breaks = 30, plot = FALSE)
hists$rand <- hist(error$rand, breaks = hists$all$breaks, plot = FALSE)
hists$fixed <- hist(error$fixed, breaks = hists$all$breaks, plot = FALSE)
hists$opt <- hist(error$opt, breaks = hists$all$breaks, plot = FALSE)
if(frac) for(i in seq_along(hists)) hists[[i]]$counts <- hists[[i]]$counts/sum(hists[[i]]$counts)
ylim <- range(c(hists$opt$counts, hists$fixed$counts, hists$rand$counts))
xlim <- c(0, max(hists$all$breaks))

labelsize <- 2.5

pdf(imgfile, width = 10, height = 10)
par(mar = c(5, 5, 0, 0), mfrow = c(3, 1))
plot(
    hists$opt, col = optcolor, border = NA,
    axes = FALSE, xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = ""
)
axis(2, cex.axis = labelsize)
grid(lty = 1)
legend("topright", bty = "n", pch = 15, col = c(optcolor, fixedcolor, randcolor), title = dynamics,
       pt.cex = 4, cex = 1.5, legend = c("Optimized", "Degree-preserving", "Random"))
plot(
    hists$fixed, col = fixedcolor, border = NA,
    axes = FALSE, xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = ""
)
axis(2, cex.axis = labelsize)
grid(lty = 1)
plot(
    hists$rand, col = randcolor, border = NA,
    axes = FALSE, xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = ""
)
axis(1, cex.axis = labelsize)
axis(2, cex.axis = labelsize)
grid(lty = 1)
dev.off()
