#' Run a simulation to assess due to unmeasured confounding
#'
#' `simBA()` runs a simulation to compute the magnitude of the bias in a hazard ratio in the presence of unmeasured confounding, possibly when proxies are available.
#'
#' @param parameters either a data.frame containing information about the data generation used in the simulation or a string containing the path to a .csv or .xlsx file containing such information. See Details for what this should contain and [create_parameters()] to create a skeleton of this object.
#' @param iterations the number of simulation iterations. Default is 500.
#' @param size the size of each sample to be generated. Default is 1000.
#' @param treatment_prevalence the desired prevalence of treatment. Should be a number between 0 and 1.
#' @param treatment_coeff the coefficient on the treatment variable in the data-generating model for survival times.
#' @param outcome_prevalence the desired prevalence of the outcome. This is used to specify a censoring time for the data-generating model.
#' @param dist the distribution to use to generate survival times. Allowable options include `"exponential"` (default) and `"weibull"`. Abbreviations allowed.
#' @param unmeasured_conf the name of the variable in `parameters` corresponding to the unmeasured confounder.
#' @param n_proxies the number of proxies for the unmeasured confounder to include in the simulation. Default is 0.
#' @param proxy_type when `n_proxies` is greater than 0, the type of variable the proxies should be. Allowable options in cldue `"binary"` (default) and `"continuous"`. Abbreviations allowed.
#' @param corr when `n_proxies` is greater than 0, the desired correlations between the proxy variable and the unmeasured confounder in the simulation. Should be length 1 (in which case all proxies have the same correlation with the unmeasured confounder) or length equal to `n_proxies`.
#' @param adj string; the method used to adjust for the confounders. Allowable options include `"matching"` (the default), which uses [MatchIt::matchit()], and `"weighting"`, which uses [WeightIt::weightit()]. Abbreviations allowed.
#' @param estimand string; the desired estimand to target. Allowable options include `"ATT"` (default), `"ATC"`, and `"ATE"`. Note this is also passed to the `estimand` argument of the function used for adjustment as specified by `adj` if omitted in `adj_args`.
#' @param adj_args a list of arguments passed to [MatchIt::matchit()] or [WeightIt::weightit()] depending on the argument to `adj`. If not supplied, the parameter defaults will be used. Take care to specify these arguments to ensure the adjustment method is as desired.
#' @param keep_data `logical`; whether to keep the datasets generated in each simulation. Default is `FALSE`. Setting to `TRUE` will make the output object large.
#' @param cl a cluster object created by [parallel::makeCluster()], or an integer to indicate number of child-processes (integer values are ignored on Windows) for parallel evaluations. See [pbapply::pbapply()] for details. Default is `NULL` for no parallel evaluation.
#' @param verbose whether to print information about the progress of the simulation, including a progress bar. Default is `TRUE`.
#'
#' @returns
#'
#' A `simBA` object, which contains the simulation outputs and their summaries. This includes the following components:
#' \describe{
#'   \item{`sim_out`}{the complete simulation results, a list with an entry for each iteration including the table of log hazard ratios, the table of standardized mean differences, and the generated dataset (if `keep_data = TRUE`)}
#'   \item{`parameters`}{the table of parameters supplied to `parameters` after some processing}
#'   \item{`SMD_table`}{the table of average standardized mean differences for the unmeasured confounder and proxies before and after matching across all iterations}
#'   \item{`HR_table`}{the table of estimated and true hazard ratios averaged across all iterations (note that log hazard ratios are averaged before exponentiating the average)}
#' }
#'
#' Basic `print()` and `summary()` methods are available. See `[plot.simBA()]` for plotting.
#'
#' @details
#'
#' `simBA()` runs a simulation study to examine the impact of an unmeasured confounder on the bias of the marginal hazard ratio when using matching or weighting to adjust for observed confounders and, optionally, proxies of the unmeasured confounder. The user must specify the simulation data-generating model using the `parameters` argument and other arguments that control generation of the treatment, outcome, and proxies. Requirements for the `parameters` input are described below. In addition, the user must specify the form of adjustment used (matching or weighting) using the `adj` argument, the desired estimand using the `estimand` argument, and any other arguments passed to `adj_args` to control the matching/weighting method. Note by default, the ATT is targeted, even though the usual default estimand for weighting using `WeightIt::weightit()` is the ATE.
#'
#' Broadly, the `parameters` input contains the name of the measured and unmeasured confounders, their variable types (binary, continuous, or count), their distributions, and their coefficients in the treatment and outcome models. These values are used to generate a synthetic dataset of size corresponding to the `size` argument, which additionally contains the true propensity score used to simulate the treatment, the treatment itself, and the outcome (i.e., survival time and whether an event occurred). When proxies are requested (i.e., `n_proxies` set to 1 or greater), proxies for the unmeasured confounder are additionally generated and appended to the synthetic dataset.
#'
#' In each iteration, a synthetic dataset is generated, and then that dataset is analyzed. First, a crude marginal hazard ratio is estimated by fitting a Cox proportional hazards model for the survival times and events as a function just of the treatment. Then, the dataset is adjusted using matching or weighting with the measured covariates, and a second hazard ratio is estimated as above, this time in the matched or weighted sample. If proxies are requested, the dataset is adjusted again using matching or weighting with the measured covariates and proxies, and a third hazard ratio is estimated as above. In addition, the balance (as measured by the  standardized mean difference \[SMD\]) is reported for the unmeasured confounder and proxies before adjustment and after each round of matching or weighting.
#'
#' ## The data-generating model
#'
#' The data-generating model for the outcome corresponds to a Cox proportional hazards model as described by Bender et al. (2005). The coefficients on the measured and unmeasured confounders in the outcome model are specified in the `parameters` input, and the coefficient on the treatment variable is specified by the `treatment_coeff` argument. The treatment is generated as a Bernoulli variable with probability equal to the true propensity score, which is generated according to a logistic regression model with the coefficients on the confounders specified in the `parameters` input.
#'
#' The proxies, if requested, are generated such that their correlation with the unmeasured confounder is exactly equal to the values supplied to `corr`. The confounder are generated as uncorrelated variables according to the distribution supplied in the `parameters` input. Binary variables are generated as Bernoulli variables with probability equal to the supplied prevalence. Continuous variables are generated as Gaussian (Normal) variables with mean and standard deviation equal to their supplied values. Count variables are generated as Poisson variables with mean equal to its supplied value.
#'
#' Some parameters are determined first by generating a dataset with one million observations. With this dataset, the intercept of the true propensity score model is selected as that which yields a treatment prevalence equal to that specified in the `treatment_prevalence` argument, and the censoring time for the outcomes is selected as that which yields an outcome event prevalence equal to that specified in the `outcome_prevalence` argument. In addition, the true marginal hazard ratio is computed using this dataset by generating potential outcomes under each treatment and fitting a Cox model of the potential outcome survival times and events as a function of the treatment under which the potential outcome was generated as recommended by Austin (2013).
#'
#' ## The `parameters` input object
#'
#' The `parameters` input must be of a specified form in order to be processed correctly. It must be a data.frame with one row for each confounder to be generated with (at least) the following columns (which are case-sensitive):
#' \describe{
#'   \item{`Variable`}{the name of the variable}
#'   \item{`Type`}{the variable type; either binary, continuous, or count (see above for how these correspond to the distribution used to generate the variable)}
#'   \item{`prevalence`}{the prevalence for binary variables (should be blank for all other variable types)}
#'   \item{`mean`}{the mean for continuous and count variable (should be blank for binary variables)}
#'   \item{`sd`}{the standard deviation for continuous variables (should be blank for all other variable types)}
#'   \item{`coeff_treatment_model`}{the coefficient on that variable in the true propensity score model for the treatment (can be blank for any variable that doens't affect treatment)}
#'   \item{`coeff_outcome_model`}{the coefficient on that variable in the outcome model for the treatment (can be blank for any variable that doens't affect the outcome)}
#' }
#'
#' The variable name supplied to `unmeasured_conf` must be present in the `parameters` input, and it must have nonzero values in both the `coeff_treatment_model` and `coeff_outcome_model` columns (or else it would not be a true confounder).
#'
#' To automatically create a skeleton of the `parameters` input for you to fill in yourself, use [create_parameters()].
#'
#'
#' @references
#' Austin PC. The performance of different propensity score methods for estimating marginal hazard ratios. *Statistics in Medicine*. 2013;32(16):2837-2849. \doi{10.1002/sim.5705}
#'
#' Bender R, Augustin T, Blettner M. Generating survival times to simulate Cox proportional hazards models. *Statistics in Medicine*. 2005;24(11):1713-1723. \doi{10.1002/sim.2059}
#'
#'
#' @seealso [create_parameters()] for creating the `parameters` input; [plot.simBA()] for plotting the results. [MatchIt::matchit()] and [WeightIt::weightit()] for the functions used for matching and weighting, respectively, which detail the defaults used by these methods and allowable arguments that can be passed to `adj_args`.
#'
#' @examples
#'
#' # Get parameters example; can also create
#' # with `create_parameters()`
#' parameters <- read.csv(system.file("extdata", "parameters.csv",
#'                                    package = "sim.BA"))
#' \donttest{
#'   # Run simulation; adjustment via PS weighting for
#'   # the ATE
#'   sim <- simBA(parameters,
#'                iterations = 50,
#'                size = 200,
#'                treatment_prevalence = .2,
#'                treatment_coef = -.25,
#'                outcome_prevalence = .5,
#'                unmeasured_conf = "u1",
#'                n_proxies = 2,
#'                proxy_type = "binary",
#'                corr = c(.5, .8),
#'                verbose = FALSE,
#'                # Adjustment arguments
#'                adj = "weighting",
#'                estimand = "ATE",
#'                adj_args = list(method = "glm"))
#'
#'   sim
#'
#'   summary(sim)
#'
#'   plot(sim, "balance")
#'
#'   plot(sim, "hr")
#' }

#' @export
simBA <- function(parameters, iterations = 500, size = 1000, treatment_prevalence,
                  treatment_coeff, outcome_prevalence, dist = "exponential",
                  unmeasured_conf, n_proxies = 0, proxy_type = "binary", corr = NULL,
                  adj = "matching", estimand = "ATT", adj_args = list(),
                  keep_data = FALSE, cl = NULL, verbose = TRUE) {

  mcall <- match.call()

  # Check arguments
  chk::chk_not_missing(unmeasured_conf, "`unmeasured_conf`")
  chk::chk_string(unmeasured_conf)

  chk::chk_count(n_proxies)

  parameters <- .process_parameters(parameters, unmeasured_conf, n_proxies)

  chk::chk_count(iterations)
  chk::chk_gt(iterations, 0)

  chk::chk_count(size)
  chk::chk_gte(size, 20)

  chk::chk_not_missing(treatment_prevalence, "`treatment_prevalence`")
  chk::chk_number(treatment_prevalence)
  chk::chk_gt(treatment_prevalence, 0)
  chk::chk_lt(treatment_prevalence, 1)

  chk::chk_not_missing(treatment_coeff, "`treatment_coeff`")
  chk::chk_number(treatment_coeff)

  chk::chk_not_missing(outcome_prevalence, "`outcome_prevalence`")
  chk::chk_number(outcome_prevalence)
  chk::chk_gt(outcome_prevalence, 0)
  chk::chk_lt(outcome_prevalence, 1)

  chk::chk_string(dist)
  dist <- tolower(dist)
  dist <- match_arg(dist, c("exponential", "weibull"))

  if (n_proxies > 0) {
    chk::chk_string(proxy_type)
    proxy_type <- tolower(proxy_type)
    proxy_type <- match_arg(proxy_type, c("binary", "continuous"))

    chk::chk_numeric(corr)
    chk::chk_gte(corr, 0)
    chk::chk_lte(corr, 1)

    if (length(corr) == 1) {
      corr <- rep(corr, n_proxies)
    }
    else if (length(corr) != n_proxies) {
      chk::err("`corr` must have length equal to `n_proxies`")
    }
  }

  chk::chk_string(adj)
  adj <- tolower(adj)
  adj <- match_arg(adj, c("matching", "weighting"))
  rlang::check_installed(switch(adj,
                                "matching" = "MatchIt",
                                "weighting" = "WeightIt"))

  chk::chk_string(estimand)
  estimand <- toupper(estimand)
  estimand <- match_arg(estimand, c("ATT", "ATC", "ATE"))

  if (length(adj_args) == 0) adj_args <- list()
  chk::chk_list(adj_args)

  if (is.null(adj_args$estimand)) {
    adj_args$estimand <- estimand
  }
  else if (!identical(adj_args$estimand, estimand)) {
    chk::wrn("the `estimand` supplied in `adj_args` does not match the estimand supplied to `estimand`")
  }

  chk::chk_flag(keep_data)

  chk::chk_flag(verbose)

  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))

  # Get params by simulating on 1e6 units
  if (verbose) {
    cat("Computing simulation parameters...\t")
  }

  sim_params <- .gen_params(parameters, size = 1e6, treatment_prevalence,
                            treatment_coeff, outcome_prevalence, dist,
                            estimand = estimand)

  if (verbose) {
    cat("Done.\nGenerating and analyzing simulated data...\n")
  }

  # Do the simulation
  sim_out <- pbapply::pblapply(seq_len(iterations), function(i) {
    simdata <- .gen_data(parameters, size, treatment_coeff, dist,
                         ps_intercept = sim_params$ps_intercept,
                         cens_time = sim_params$cens_time,
                         true_logHR = sim_params$true_logHR,
                         unmeasured_conf, n_proxies, proxy_type, corr)

    .analyze_data(simdata, keep_data = keep_data,
                  adj = adj, adj_args = adj_args)
  })

  balance_table <- do.call("rbind", lapply(sim_out, `[[`, "all_balance_tables"))

  results_n <- do.call("rbind", lapply(sim_out, `[[`, "all_HRs"))

  ## Prepare the the SMD table to output
  balance_table$SMD[is.na(balance_table$SMD)] <- 0

  SMD_table <- aggregate(SMD ~ variable + adjustment,
                         data = balance_table, FUN = mean)

  SMD_table <- SMD_table[order(SMD_table$variable, SMD_table$adjustment),]

  ## Prepare the HR table to output
  HR_table <- aggregate(HR ~ adjustment,
                         data = results_n, FUN = mean)

  HR_table <- HR_table[order(HR_table$adjustment),]

  HR_table$HR <- exp(HR_table$HR)

  out <- list(sim_out = sim_out,
              parameters = parameters,
              SMD_table = SMD_table,
              HR_table = HR_table)

  if (verbose) {
    cat("Done.\n")
  }

  attr(out, "size") <- size
  attr(out, "treatment_prevalence") <- treatment_prevalence
  attr(out, "treatment_coeff") <- treatment_coeff
  attr(out, "outcome_prevalence") <- outcome_prevalence
  attr(out, "dist") <- dist
  attr(out, "unmeasured_conf") <- unmeasured_conf
  attr(out, "proxy_type") <- proxy_type
  attr(out, "corr") <- corr
  attr(out, "adj") <- adj
  attr(out, "call") <- mcall

  class(out) <- "simBA"

  out
}

#' @exportS3Method print simBA
print.simBA <- function(x, ...) {
  cat("A `simBA` object\n")
  cat(sprintf(" - sample size: %s\n", attr(x, "size")))
  cat(sprintf(" - treatment prevalence: %s\n", attr(x, "treatment_prevalence")))
  cat(sprintf(" - treatment coefficient: %s (true HR: %s)\n",
              attr(x, "treatment_coeff"),
              round(x$HR_table$HR[x$HR_table$adjustment == "true"], 2)))
  cat(sprintf(" - outcome prevalence: %s\n", attr(x, "outcome_prevalence")))
  cat(sprintf(" - outcome distribution: %s\n", firstup(attr(x, "dist"))))
  cat(sprintf(" - unmeasured confounder: %s\n", add_quotes(attr(x, "unmeasured_conf"), 2)))
  n_proxies <- length(attr(x, "corr"))
  if (n_proxies == 0) {
    cat(" - proxies: none\n")
  }
  else {
    cat(sprintf(" - proxies: %s (%s; r = %s)\n",
                n_proxies,
                attr(x, "proxy_type"),
                paste(round(unique(attr(x, "corr")), 2), collapse = ", ")))
  }
  cat(sprintf(" - adjustment method: %s\n", attr(x, "adj")))
  cat("Use `summary()` and `plot()` to examine results.\n")

  invisible(x)
}

#' @exportS3Method summary simBA
summary.simBA <- function(object, ...) {
  out <- object[c("SMD_table", "HR_table")]

  class(out) <- "summary.simBA"

  out
}

#' @exportS3Method print summary.simBA
print.summary.simBA <- function(x, digits = 3, ...) {
  names(x$SMD_table) <- firstup(names(x$SMD_table))
  cat("- Average Standardized Mean Difference (SMD):\n")
  print.data.frame(x$SMD_table, row.names = FALSE, digits = digits)

  names(x$HR_table) <- firstup(names(x$HR_table))
  cat("\n- Hazard Ratios (HR):\n")
  print.data.frame(x$HR_table, row.names = FALSE, digits = digits)

  invisible(x)
}

.gen_params <- function(parameters, size = 1e6, treatment_prevalence, treatment_coeff,
                        outcome_prevalence, dist = "exponential",
                        estimand = "ATT") {

  # Generate PS intercept, censoring time, and true marginal HR

  size <- switch(estimand,
                 "ATT" = ceiling(size/treatment_prevalence),
                 "ATC" = ceiling(size/(1 - treatment_prevalence)),
                 size)

  ## Generate simuldation data matrix
  simdata <- vapply(seq_len(nrow(parameters)), function(i) {
    switch(parameters$Type[i],
           "binary" = rbinom(size, 1, parameters$prevalence[i]),
           "continuous" = rnorm(size, parameters$mean[i], parameters$sd[i]),
           "count" = rpois(size, parameters$mean[i])
    )}, numeric(size))

  colnames(simdata) <- parameters$Variable

  ## Generate linear predictor of the true PS model for treatment
  treat_linear_pred <- drop(simdata %*% parameters$coeff_treatment_model)

  ## Compute PS model intercept to match requested treatment prevalence

  prob.exposure <- runif(size)

  optfun <- function(intercept) {
    trueps <- plogis(intercept + treat_linear_pred)
    treatment <- as.integer(trueps > prob.exposure)
    generated_treatment_prev <- sum(treatment) / size

    # Minimize square of treatment prevalence discrepancy
    (generated_treatment_prev - treatment_prevalence) ^ 2
  }

  opt_out <- optimize(optfun, interval = c(-40, 40))
  intercept <- opt_out$minimum

  trueps <- plogis(intercept + treat_linear_pred)
  treatment <- as.integer(trueps > prob.exposure)

  ## Generate linear predictor of the outcome; treatment effect
  ## will be added later
  surv_linear_pred <- drop(simdata %*% parameters$coeff_outcome_model)

  Uevent <- runif(size)

  scale.event <- 0.0001 # scale parameter for exponentially distributed event time
  shape.event <- switch(dist, "exponential" = 1, "weibull" = 2)

  ## Building the outcome time-to-event model
  .gen_true_time <- function(linear_pred, Uevent) {
    (-log(Uevent)/(scale.event * exp(linear_pred)))^(1/shape.event)
  }

  ## Draw an event time and a censoring time, the minimum  of these is "observed"
  ## and we record whether it was the event time or the censoring time

  true_time <- .gen_true_time(surv_linear_pred + treatment * treatment_coeff, Uevent)

  ## Censor at nth percentile of the generated event time to generate desired
  ## outcome prevalence
  cens_time <- quantile(true_time, outcome_prevalence)

  ## Generate true potential outcomes to estimate true marginal log HR
  true_time0 <- .gen_true_time(surv_linear_pred, Uevent)
  true_time1 <- .gen_true_time(surv_linear_pred + treatment_coeff, Uevent)

  dat <- data.frame(.true_time = c(true_time0, true_time1),
                    .treatment = rep(0:1, each = size),
                    .true_treatment = c(treatment, treatment))

  dat$.true_survt <- pmin(cens_time, dat$.true_time)
  dat$.event <- dat$.true_survt == dat$.true_time

  ## Subset data to get POs for given estimand
  dat <- switch(estimand,
                "ATT" = dat[dat$.true_treatment == 1,],
                "ATC" = dat[dat$.true_treatment == 0,],
                dat)

  ## Fit Cox model on POs and synthetic treatments; seee Austin (2013)
  fit <- survival::coxph(survival::Surv(.true_survt, .event) ~ .treatment,
                         data = dat,
                         robust = FALSE)

  ## Return PS intercept, censoting time, and true marginal logHR
  list(ps_intercept = intercept,
       cens_time = cens_time,
       true_logHR = fit$coefficients[".treatment"])
}

.gen_data <- function(parameters, size,
                      treatment_coeff,
                      dist = "exponential",
                      ps_intercept, cens_time,
                      true_logHR,
                      unmeasured_conf,
                      n_proxies = 0, proxy_type, corr = NULL) {

  # Generate covariates with pre-specified distribution
  simdata <- vapply(seq_len(nrow(parameters)), function(i) {
    switch(parameters$Type[i],
           "binary" = rbinom(size, 1, parameters$prevalence[i]),
           "continuous" = rnorm(size, parameters$mean[i], parameters$sd[i]),
           "count" = rpois(size, parameters$mean[i])
    )}, numeric(size))

  colnames(simdata) <- parameters$Variable

  covariates <- setdiff(parameters$Variable, unmeasured_conf)

  # Generate treatment and outcome

  ## Treatment model linear predictor
  treat_linear_pred <- drop(simdata %*% parameters$coeff_treatment_model)

  ## Add intercept and take plogis() to get true PS
  trueps <- plogis(ps_intercept + treat_linear_pred)

  ## Generate treatment as binomial variable
  treatment <- rbinom(size, 1, trueps)

  ## Generate linear predictor of the outcome; treatment effect
  ## will be added later
  surv_linear_pred <- drop(simdata %*% parameters$coeff_outcome_model)

  Uevent <- runif(size)
  scale.event <- 0.0001 # scale parameter for exponentially distributed event time
  shape.event <- switch(dist, "exponential" = 1, "weibull" = 2)

  ## Draw an event time and a censoring time, the minimum  of these is "observed"
  ##  and we record whether it was the event time or the censoring time
  .gen_true_time <- function(linear_pred, Uevent) {
    (-log(Uevent)/(scale.event * exp(linear_pred)))^(1/shape.event)
  }

  true_time <- .gen_true_time(surv_linear_pred + treatment * treatment_coeff, Uevent)

  ## Observed time is min of censored and true
  survt <- pmin(cens_time, true_time)

  ## Set event to 1 if event is observed
  event <- survt == true_time

  ## Bring together in data.frame
  simdata <- as.data.frame(simdata)
  simdata$.trueps <- trueps
  simdata$.treatment <- treatment
  simdata$.survt <- survt
  simdata$.event <- event

  # Insert proxies if requested
  if (n_proxies > 0) {

    ## Add proxies withspecified correlation
    add_proxy <- {
      if (proxy_type == "binary") {
        function(x, rho) {

          #Calculate the mean of x
          x_bar <- mean(x)

          #Calculate the expected value of xy under the assumption that x and y are correlated with correlation coefficient rho
          xy_bar <- rho * x_bar + (1 - rho) * x_bar ^ 2

          #Calculate the number of 1s in x that need to be flipped to achieve the desired correlation coefficient
          toflip <- sum(x == 1) - round(length(x) * xy_bar)

          #Create a copy of x
          y <- x

          #Randomly flip the values of x at the specified indexes to achieve the desired correlation coefficient
          y[sample(which(x == 0), toflip)] <- 1
          y[sample(which(x == 1), toflip)] <- 0

          y
        }
      }
      else {
        function(x, rho) {
          (rho * (x - mean(x)))/sqrt(var(x)) + sqrt(1 - rho^2) * rnorm(length(x))
        }
      }
    }

    proxy_vars <- paste0("p", seq_len(n_proxies))

    for (i in seq_len(n_proxies)) {
      simdata[[proxy_vars[i]]] <- add_proxy(simdata[[unmeasured_conf]], corr[i])
    }
  }
  else {
    proxy_vars <- character(0)
  }

  attr(simdata, "covariates") <- covariates
  attr(simdata, "proxy_vars") <- proxy_vars
  attr(simdata, "unmeasured_conf") <- unmeasured_conf
  attr(simdata, "true_logHR") <- true_logHR

  simdata
}

.analyze_data <- function(simdata, keep_data = FALSE, adj = "matching", adj_args = list()) {
  # Note: robust SEs set to FALSE, as they are not used in the simulation. Correct use
  # requires including subclass membership for matching, too.

  adj_fun <- {
    if (adj == "matching")
      MatchIt::matchit
    else if (length(adj_args) > 1 || !adj_args$estimand %in% c("ATT", "ATC", "ATE", "ATO"))
      WeightIt::weightit
    else
      #Simple PS weighting with no dependencies
      function(formula, data, estimand) {
        fit <- glm(formula, data = data, family = quasibinomial)
        ps <- fit$fitted.values
        w_ATE <- ifelse(simdata$.treatment == 1, 1 / ps, 1 / (1 - ps))

        list(weights = switch(
          estimand,
          "ATT" = ps * w_ATE,
          "ATC" = (1 - ps) * w_ATE,
          "ATO" = ps * (1 - ps) * w_ATE,
          "ATE" = w_ATE)
        )
      }
  }

  adj_args$data <- simdata


  ## true effects
  true_results <- data.frame(HR = attr(simdata, "true_logHR"),
                             SE = NA_real_,
                             adjustment = "true")

  ## model 1- crude
  crude_model <- survival::coxph(survival::Surv(.survt, .event) ~ .treatment,
                                 data = simdata,
                                 robust = FALSE,
                                 method = "breslow")
  summ <- summary(crude_model)$coefficients

  HR <- summ[".treatment", "coef"]
  SE <- NA #summ[".treatment", "robust se"]

  crude_results <- data.frame(HR, SE, adjustment = "crude")

  ## model 2- adjusted for all measured confounding (L1)
  covariates <- attr(simdata, "covariates")

  f_L1 <- reformulate(covariates, ".treatment")

  adj_args$formula <- f_L1

  adj.out_L1 <- do.call(adj_fun, adj_args, quote = TRUE)

  simdata$.L1_weights <- adj.out_L1$weights
  simdata$.L1_subclass <- adj.out_L1$subclass

  L1_model <- survival::coxph(survival::Surv(.survt, .event) ~ .treatment,
                              data = simdata,
                              weights = simdata$.L1_weights,
                              subset = simdata$.L1_weights != 0,
                              # cluster = matched_L2$subclass,
                              robust = FALSE,
                              method = "breslow")

  summ <- summary(L1_model)$coefficients

  HR <- summ[".treatment", "coef"]
  SE <- NA #summ[".treatment", "robust se"]
  L1_results <- data.frame(HR, SE, adjustment = "L1")

  all_HRs <- rbind(true_results, crude_results, L1_results)

  ## model 3- adjusted for all measured confounding (L1) plus the proxies
  proxy_vars <- attr(simdata, "proxy_vars")

  if (length(proxy_vars) > 0) {

    # write the formula
    f_L2 <- reformulate(c(covariates, proxy_vars), ".treatment")

    adj_args$formula <- f_L2

    adj.out_L2 <- do.call(adj_fun, adj_args, quote = TRUE)

    simdata$.L2_weights <- adj.out_L2$weights
    simdata$.L2_subclass <- adj.out_L2$subclass

    L2_model <- survival::coxph(survival::Surv(.survt, .event) ~ .treatment,
                                data = simdata,
                                weights = simdata$.L2_weights,
                                subset = simdata$.L2_weights != 0,
                                # cluster = matched_L2$subclass,
                                robust = FALSE,
                                method = "breslow")

    summ <- summary(L2_model)$coefficients

    HR <- summ[".treatment", "coef"]
    SE <- NA #summ[".treatment", "robust se"]
    L2_results <- data.frame(HR, SE, adjustment = "L2")

    all_HRs <- rbind(all_HRs, L2_results)
  }

  all_HRs$adjustment <- factor(all_HRs$adjustment, levels = c("crude", "L1", "L2", "true"))

  #Balance
  unmeasured_conf <- attr(simdata, "unmeasured_conf")

  if (length(proxy_vars) == 0) {
    balance_vars <- unmeasured_conf

    #Crude
    balance_table_crude <- data.frame(variable = balance_vars,
                                      adjustment = "crude",
                                      SMD = cobalt::col_w_smd(simdata[balance_vars],
                                                              treat = simdata$.treatment,
                                                              s.d.denom = "treated"))

    #L1
    balance_table_L1 <- data.frame(variable = balance_vars,
                                   adjustment = "L1",
                                   SMD = cobalt::col_w_smd(simdata[balance_vars],
                                                           treat = simdata$.treatment,
                                                           s.d.denom = "treated",
                                                           weights = simdata$.L1_weights))

    all_balance_tables <- rbind(balance_table_crude, balance_table_L1)
  }
  else {
    balance_vars <- c(unmeasured_conf, proxy_vars)

    #Crude
    balance_table_crude <- data.frame(variable = balance_vars,
                                      adjustment = "crude",
                                      SMD = cobalt::col_w_smd(simdata[balance_vars],
                                                              treat = simdata$.treatment,
                                                              s.d.denom = "treated"))

    #L1
    balance_table_L1 <- data.frame(variable = balance_vars,
                                   adjustment = "L1",
                                   SMD = cobalt::col_w_smd(simdata[balance_vars],
                                                           treat = simdata$.treatment,
                                                           s.d.denom = "treated",
                                                           weights = simdata$.L1_weights))

    #L2
    balance_table_L2 <- data.frame(variable = balance_vars,
                                   adjustment = "L2",
                                   SMD = cobalt::col_w_smd(simdata[balance_vars],
                                                           treat = simdata$.treatment,
                                                           s.d.denom = "treated",
                                                           weights = simdata$.L2_weights))

    all_balance_tables <- rbind(balance_table_crude, balance_table_L1, balance_table_L2)
  }

  all_balance_tables$adjustment <- factor(all_balance_tables$adjustment, levels = c("crude", "L1", "L2"))

  out <- list(all_HRs = all_HRs,
              all_balance_tables = all_balance_tables)

  if (keep_data) {
    out$data <- simdata
  }

  out
}
