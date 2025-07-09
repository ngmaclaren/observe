library(parallel)
library(nlme)

get_testerrors <- function(datalist) {
    errors <- datalist$testerror
    sapply(errors, mean)
}

get_error_df <- function(dynamics, TAPAS = FALSE) {
    if(TAPAS) {
        dl <- readRDS(paste0("../data/datalist-TAPAS-fmri-", dynamics, ".rds"))
        flag <- readRDS("./TAPAS-flag4.rds")
        ##flag <- sapply(dl, `[[`, "flag")
    } else {
        dl <- readRDS(paste0("../data/datalist-fmri-", dynamics, ".rds"))
        flag <- sapply(dl, `[[`, "flag")
    }
    dlr <- dl[which(flag == "use")]
    error <- data.frame(t(sapply(dlr, get_testerrors)))
    error$dynamics <- dynamics
    error$pID <- seq_len(nrow(error))
    error
}

dyns <- c("doublewell", "mutualistic", "SIS", "genereg", "wilsoncowan") # , "SIS"
dfs <- mclapply(dyns, get_error_df, TAPAS = TRUE, mc.cores = 5)
df <- do.call(rbind, dfs)

rdf <- reshape(
    df,
    varying = c("opt", "fixed", "rand"),
    v.names = "error",
    timevar = "ns.type",
    times = c("opt", "fixed", "rand"),
    new.row.names = seq_len(nrow(df)*3),
    direction = "long"
)
rdf$dynamics <- factor(rdf$dynamics, levels = dyns)
rdf$ns.type <- factor(rdf$ns.type, levels = c("rand", "opt", "fixed"))

model <- lme(error ~ dynamics + ns.type, random = ~ 1 | pID, data = rdf)
summary(model)

with(list(x = summary(model)$tTable), x["ns.typeopt", "Value"]/x["(Intercept)", "Value"])
with(list(x = summary(model)$tTable), x["ns.typefixed", "Value"]/x["(Intercept)", "Value"])
## with(list(x = summary(model)$tTable), x["dynamicsgenereg", "Value"]/x["(Intercept)", "Value"])
## with(list(x = summary(model)$tTable), x["dynamicswilsoncowan", "Value"]/x["(Intercept)", "Value"])

## notapprop <- subset(rdf, dynamics %in% c("doublewell", "mutualistic", "SIS"))

## notapprop_model <- lme(error ~ dynamics + ns.type, random = ~ 1 | pID, data = notapprop)
## summary(notapprop_model)

## with(list(x = summary(notapprop_model)$tTable), x["ns.typeopt", "Value"]/x["(Intercept)", "Value"])
## with(list(x = summary(notapprop_model)$tTable), x["ns.typefixed", "Value"]/x["(Intercept)", "Value"])

## approp <- subset(rdf, dynamics == "genereg" | dynamics == "wilsoncowan")

## approp_model <- lme(error ~ dynamics + ns.type, random = ~ 1 | pID, data = approp)
## summary(approp_model)

## with(list(x = summary(approp_model)$tTable), x["ns.typeopt", "Value"]/x["(Intercept)", "Value"])
## with(list(x = summary(approp_model)$tTable), x["ns.typefixed", "Value"]/x["(Intercept)", "Value"])

model_dyn <- function(dyn) {
    sdf <- subset(rdf, dynamics == dyn)
    return(lme(error ~ ns.type, random = ~ 1 | pID, data = sdf))
}

dynmodels <- lapply(dyns, model_dyn)
names(dynmodels) <- dyns
lapply(dynmodels, summary)

for(i in seq_along(dynmodels)) {
    print(names(dynmodels)[i])
    m <- dynmodels[[i]]
    print(with(list(x = summary(m)$tTable),
               paste0(round(100*x["ns.typeopt", "Value"]/x["(Intercept)", "Value"], 2), "%")))
    print(with(list(x = summary(m)$tTable),
               paste0(round(100*x["ns.typefixed", "Value"]/x["(Intercept)", "Value"], 2), "%")))
}; rm(m)

lapply(dynmodels, function(x) summary(x)$tTable)

## remove_outliers <- function(sdf) {
##     outliers <- unique(unlist(
##         apply(sdf[, 1:3], 2, function(x) which(x >= mean(x) + 4*sd(x)), simplify = FALSE)
##     ))
##     return(sdf[-outliers, ])
## }

## sdfs <- split(df, df$dynamics)
## sdfs <- lapply(sdfs, remove_outliers)
## df <- do.call(rbind, sdfs)

