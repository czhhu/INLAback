#'@name INLAback
#'@param fam1 String defining the likelihood familiy.
#'@param dataf A dataframe including covariates and response data.
#'@param invariant The part of the formula that should not change
#'@param direction string 'backwards' for backwards variable elimination.
#'@param include Vector of integers to determine which columns in dataf should be used. If NULL, use all columns except y.
#'@param y String determining the response variable.
#'@param ... Further arguments to \code{INLA::inla} function.
#'@importFrom stats formula
#'@export



INLAback<-function(fam1 = "gaussian",
                   dataf,
                   invariant = "Intercept",
                   direction = "backwards",
                   y = NULL,
                   include = NULL,
                   powerl = 1,
                   inter = 1,
                   thresh = 2,
                   num.threads = 1,
                   ...) {



  # Basic checks
  if (is.null(nrow(dataf))) {
    stop("error in data frame")
  }
  if (nrow(dataf) == 0) {
    stop("no rows in data frame")
  }
  if (is.null(y)) {
    stop("no y variable")
  }
  if (class(dataf) == "data.frame") {
    stop("data is not a data frame")
  }



  if (is.null(include))
    include <- (1:ncol(dataf))[!names(dataf) %in% y]


  z <- NULL


  facts <- sapply(dataf, is.factor)[include]
  explF<-names(dataf)[include]
  expl<-explF[!facts]


  expl <- expandExplanatoryVars(expl, explF, facts)

  choice <- NULL
  chosen <- NULL
  new1 <- NULL
  dicloss <- 999
  dicold <- NULL


  while(length(expl) > 0){
    if (direction == "backwards") {
      runs <- c(1:length(expl), 9999999)
    } else{
      runs <- 1:length(expl)
    }

    for (ii in runs) {
      if (direction == "backwards") {
        if (ii == 9999999) {
          ii <- 1:length(expl)
        } else{
          ii <-
            {
              -1 * ii
            } # Drop each variable in turn. Final run is ii == 9999 and therefore do all variables.
        }
      }

      if (is.null(chosen)) {
        if (length(expl[ii]) > 0) {
          formula2 <-
            formula(paste(y, "~", invariant, "+", paste(expl[ii], collapse = "+"), sep =
                            ""))
        } else {
          formula2 <- formula(paste(y, "~", invariant))
        }
      } else{
        if (length(expl[ii]) > 0) {
          formula2 <-
            formula(paste(y, "~", invariant, "+", chosen, " + ", expl[ii], sep = ""))
        } else {
          formula2 <- formula(paste(y, "~", invariant, "+", chosen))
        }
      }

      result2 <- INLA::inla(
        formula2,
        family = fam1,
        num.threads = num.threads,
        control.compute = list(cpo = TRUE,
                               dic = TRUE,
                               waic = TRUE),
        verbose = FALSE,
        data = dataf,
        ...
      )

      rmse <-
        sqrt(mean((
          dataf[, y] - result2$summary.fitted.values$mean[1:nrow(dataf)]
        ) ^ 2, na.rm = TRUE))
      sumcpo <- sum(log(result2$cpo$cpo), na.rm = TRUE)
      if (length(ii) > 1) {
        var1 <- paste(expl[ii], collapse = "+")
      } else{
        var1 <- expl[abs(ii)]
      }
      if (is.null(choice)) {
        choice <-
          data.frame(
            var = var1,
            aic = result2$waic$waic,
            rmse,
            sumcpo,
            stringsAsFactors = FALSE
          )
      } else{
        choice <-
          rbind(
            choice,
            data.frame(
              var = var1,
              aic = result2$waic$waic,
              rmse,
              sumcpo,
              stringsAsFactors = FALSE
            )
          )
      }
    }##end of run through

    new1 <- choice[which.min(choice$aic), 1]
    # If not the first time through, calculate dic loss
    if (!is.null(dicold)) {
      dicloss <- dicold - min(choice$aic, na.rm = TRUE)[1]
    }
    # Update dic old
    dicold <- choice[which.min(choice$aic), 2]

    if (is.null(z)) {
      progress <- choice[choice$var == new1, ]
      z <- 1
    } else {
      progress <- rbind(progress, choice[choice$var == new1, ])
    }

    message(paste(new1, " - ", min(choice$aic, na.rm = TRUE)), sep = "")
    choice <- NULL
    if (dicloss > thresh) {

      if (direction == "backwards") {
        expl <- expl[!expl == new1]
      }
      if (direction == "forwards") {
        if (is.null(chosen)) {
          chosen <- new1
          expl <- expl[!expl == new1]
        } else {
          chosen <- paste(chosen, " + ", new1, sep = "")
          expl <- expl[!expl == new1]
        }
      }
    } else {
      break
    }

  }

  if (direction == "backwards") {
    formulax <-
      formula(paste(y, "~", invariant, "+", paste(expl[ii], collapse = "+"), sep = ""))
  } else {
    if(!is.null(chosen)){
      formulax <- formula(paste(y, "~", invariant, "+", chosen, sep = ""))
    } else {
      formulax <- formula(paste(y, "~", invariant, sep = ""))
    }
  }
  #print(formulax)

  result2 <- INLA::inla(
    formulax,
    family = fam1,
    num.threads = num.threads,
    control.compute = list(cpo = TRUE,
                           dic = TRUE,
                           waic = TRUE),
    verbose = FALSE,
    data = dataf,
    ...
  )


  output <- list(
    best_formula = formulax,
    waic = dicold,
    progress = progress,
    best_model = result2
  )
  class(output) <- 'INLAback'

  return(output)

}##end of function
