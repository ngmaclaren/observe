library(optparse)
optionlist <- list(
    make_option(
        "--dynamics", type = "character", default = "doublewell",
        help = "The dynamics on the network. Default is %default. Options: 'doublewell', 'SIS', 'genereg', 'mutualistic', and 'wilsoncowan'."
    ),
    make_option(
        "--test", action = "store_true", default = FALSE,
        help = "If flag is present, don't run the simulation for the whole data set. Instead, select an arbitrary individual's data and run on that data only."
    )
)
args <- parse_args(
    OptionParser(option_list = optionlist),
    convert_hyphens_to_underscores = TRUE
)

library(parallel)
ncores <- detectCores()-2
library(igraph)
library(sdn)
## library(optNS)
source("~/Documents/reduction/sims/local-optNS.R")

.wilsoncowan <- list(
    xinit.low = 0, xinit.high = 8,
    mu = 3, delta = 1,
    u = 0, use.matrix = TRUE, # because all fMRI networks are weighted
    D = 0.05, Ds = seq(0, 2, length.out = 100)
)
wilsoncowan <- function(t, x, params) {
    with(params, {
        coupling <- D*rowSums(A*outer(rep(1, length(x)), 1/(1 + exp(mu - delta*x))))
        dx <- -x + coupling + u
        return(list(c(dx)))
    })
}

if(interactive()) {
    args$dynamics <- "doublewell"
    args$test <- TRUE
}

test <- args$test
dynamics <- args$dynamics

outfile <- paste0("datalist-TAPAS-nonsparse-", dynamics, ".rds")
print(outfile)

model <- get(dynamics)
params <- get(paste0(".", dynamics))

## datadir <- "~/Documents/reduction/data/TAPAS/"
datadir <- "/projects/academic/naokimas/neil/TAPASnonsparse/"
datafiles <- list.files(path = datadir, pattern = ".csv")

tsdir <- "/projects/academic/naokimas/neil/brains-ns50/"
tsfiles <- list.files(path = tsdir, pattern = ".txt")

stopifnot(all.equal(gsub(".txt", "", tsfiles), gsub("matrix.csv", "", datafiles)))

make_nodesets <- function(datafile, tsfile) {
    df <- read.table(paste0(tsdir, tsfile), sep = " ")
    nr <- nrow(df)

    X.train <- as.matrix(df[seq_len(nr-100), ])
    X.test <- as.matrix(df[seq(nr-100+1, nr, by = 1), ])

    avgxi <- colMeans(X.test)
    mu <- mean(avgxi) 
    sigma <- sd(avgxi) 
    threshold <- mu + 5*sigma 
    test <- any(abs(avgxi) >= threshold) 
    if(isFALSE(test)) flag <- "use" else flag <- "dontuse"

                                        # The TAPAS project uses j -> i rather than i -> j
                                        # as the adjacency matrix convention.
    A <- t(as.matrix(read.csv(paste0(datadir, datafile), header = FALSE)))

                                        # These networks have self-loops
    g <- graph_from_adjacency_matrix(A, "directed", weighted = TRUE)#, diag = FALSE) 
    N <- vcount(g)
    n <- floor(log(N))

    times <- 0:15
    control <- list(times = times, ncores = ncores) # 1
    params <- c(params, list(A = A))
    params$use.matrix <- TRUE
    Ds <- seq(0, 1, length.out = 100) # this is right for the Pearson networks
    ## us <- switch(
    ##     dynamics,
    ##     mutualistic = seq(0, 1, length.out = 100),
    ##     seq(0, 5, length.out = 100)
    ## )
    
    xinit <- switch(
        dynamics,
        genereg = params$xinit.high,
        params$xinit.low
    )

    Y <- solve_in_range(
        Ds, "D", model, rep(xinit, N), params, control, "ode",
        ## us, "u", model, rep(xinit, N), params, control, "ode",
        method = "adams", maxsteps = 100000
    )
    y <- rowMeans(Y)
    print(datafile)
    print(dim(Y))
    range(Y)
    ## y
    ## apply(Y, 1, min)
    ## apply(Y, 1, function(x) sum(x < 0))

    opts <- make_dataset(
        ntrials = 100, ns.type = "opt", ncores = ncores, n = n, g = g, y = y, Y = Y
    )

    best <- opts[[which.min(get_error(opts))]]

    fixeds <- make_dataset(
        ntrials = 100, ns.type = "fixed", ncores = ncores, n = n, g = g, comps = best$vs, y = y, Y = Y
    )

    rands <- make_dataset(
        ntrials = 100, ns.type = "rand", ncores = ncores, n = n, g = g, y = y, Y = Y
    )

    testerror <- list(
        opt = sapply(opts, function(ns) obj_fn(ns$vs, rowMeans(X.test), X.test)),
        fixed = sapply(fixeds, function(ns) obj_fn(ns$vs, rowMeans(X.test), X.test)),
        rand = sapply(rands, function(ns) obj_fn(ns$vs, rowMeans(X.test), X.test))
    )

    return(
        list(
            ##network = g,
            opts = opts,
            fixeds = fixeds,
            rands = rands,
            testerror = testerror,
            flag = flag
        )
    )
}

if(test) {
    datafile <- datafiles[14]
    tsfile <- tsfiles[14]
    out <- make_nodesets(datafile, tsfile)
} else {
    datalist <- mapply(make_nodesets, datafiles, tsfiles, SIMPLIFY = FALSE)
    pIDs <- gsub("matrix.csv", "", datafiles)
    names(datalist) <- pIDs
    saveRDS(datalist, outfile)
}

### old code

    
## sapply(testerror, summary)
##A <- as.matrix(read.csv(paste0("../data/", pID, "matrix.csv"), header = FALSE))
## pID <- "100307" # 100206

## datadir <- "/projects/academic/naokimas/neil/brains-ns50/"
## datafile <- paste0(pID, ".txt")
## dynamics <- "doublewell"
## model <- get(dynamics)
## params <- get(paste0(".", dynamics))
