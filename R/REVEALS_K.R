library(disqover)

REVEALSinR_K <- function (pollen, params, tBasin, dBasin, dwm = "lsm unstable", 
          n = 1000, regionCutoff = 1e+05, ppefun = rnorm_reveals, pollenfun = rmultinom_reveals, 
          writeresults = FALSE, verbose = TRUE) 
{
  dwm <- tolower(dwm)
  tBasin <- tolower(tBasin)
  names(params) <- tolower(names(params))
  pollen <- as.data.frame(pollen)
  params <- as.data.frame(params)
  row.names(params) <- params[, 1]
  params <- params[, -1]
  if (!checkREVEALS(pollen = pollen, params = params, tBasin = tBasin, 
                    dBasin = dBasin, dwm = dwm, n = n, regionCutoff = regionCutoff, 
                    writeresults = writeresults, verbose = verbose)) 
    stop("Some parameter is not suitable - check and correct")
  params <- params[colnames(pollen[, -1]), ]
  ppes <- params[, c("ppe", "ppe.error")]
  vg <- params[, "fallspeed"]
  if (verbose) {
    message("REVEALS:\tCalculating deposition factor for each species ...")
  }
  deposition <- do.call(rbind, lapply(vg, DispersalFactorK, 
                                      tBasin = tBasin, dBasin = dBasin, dwm = dwm, regionCutoff = regionCutoff))
  all <- apply(pollen, 1, REVEALS, n = n, ppes = ppes, ppefun = ppefun, 
               deposition = deposition, pollenfun = pollenfun, writeresults = writeresults, 
               verbose = verbose)
  results_df <- list(data.frame(dwm, tBasin, pollen[1], t(sapply(all, 
                                                            FUN = function(x) x$meansim)), t(sapply(all, FUN = function(x) x$mediansim)), 
                           t(sapply(all, FUN = function(x) x$q90sim)), t(sapply(all, 
                                                                                FUN = function(x) x$q10sim)), t(sapply(all, FUN = function(x) x$sdsim))),
                     deposition)
  colnames(results_df[[1]]) <- c("Distance.weighting", "Basin.type", 
                            names(pollen)[1], paste(names(pollen[-1]), rep(c("mean", 
                                                                             "median", "q90", "q10", "sd"), each = length(pollen[1, 
                                                                                                                                 ]) - 1), sep = "."))
  if (writeresults) {
    nn <- names(results_df)
    k <- c(1:4, grep(x = nn, pattern = "mediansim"), grep(x = nn, 
                                                          pattern = "q10sim"), grep(x = nn, pattern = "q90sim"))
    write.table(results_df[, k], file = "ResultStatistics_confidence_limits.csv", 
                append = FALSE, row.names = FALSE, col.names = TRUE, 
                sep = ";", dec = ",")
    k <- c(1:4, grep(x = nn, pattern = "meansim"), grep(x = nn, 
                                                        pattern = "sdsim"))
    write.table(results_df[, k], file = "ResultStatistics_SD.csv", 
                append = FALSE, row.names = FALSE, col.names = TRUE, 
                sep = ";", dec = ",")
  }
  invisible(results_df)
}