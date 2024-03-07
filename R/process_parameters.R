.process_parameters <- function(parameters, unmeasured_conf, n_proxies = 0) {
  if (is.character(parameters)) {
    if (length(parameters) != 1) {
      chk::err("`parameters` should be a data.frame or the name of a file containing the parameters")
    }

    if (endsWith(parameters, ".csv")) {
      parameters <- utils::read.csv(parameters)
    }
    else if (endsWith(parameters, ".xlsx")) {
      rlang::check_installed("openxlsx")

      parameters <- as.data.frame(openxlsx::read.xlsx(parameters, 1))
    }
    else {
      chk::err("if supplied as a string, `parameters` must refer to a .csv or .xlsx file")
    }

    df_name <- "the data.frame in the file named in `parameters`"
  }
  else if (is.data.frame(parameters)) {
    #Un-tibble
    parameters <- as.data.frame(parameters)

    df_name <- "`parameters`"
  }
  else {
    chk::err("`parameters` should be a data.frame or the name of a file containing the parameters")
  }

  required_names <- c("Variable", "Type", "prevalence", "mean", "sd",
                      "coeff_treatment_model", "coeff_outcome_model")
  if (is.null(names(parameters)) || ncol(parameters) < length(required_names) ||
      !all(required_names %in% names(parameters))) {
    chk::err(sprintf("%s should have the following column names: %s. Use `create_parameters()` to create it",
                     df_name, paste(required_names, collapse = ", ")))
  }

  if (nrow(parameters) == 0) {
    chk::err(sprintf("%s should have at least one row", df_name))
  }

  #Variable
  parameters$Variable <- as.character(parameters$Variable)
  if (anyNA(parameters$Variable) || any(parameters$Variable == "")) {
    chk::err(sprintf("missing and empty values are not allowed in the `Variable` column in %s", df_name))
  }

  if (!unmeasured_conf %in% parameters$Variable) {
    chk::err(sprintf("the name supplied to `unmeasured_conf` (\"%s\") must be present in the `Variable` column in %s",
                     unmeasured_conf, df_name))
  }

  if (!chk::vld_unique(parameters$Variable)) {
    chk::err(sprintf("duplicate values are not allowed in the `Variable` column in %s",
                     df_name))
  }

  if (any(startsWith(parameters$Variable, "."))) {
    chk::err(sprintf("the names in the `Variable` column in %s cannot start with a period/full stop ('.')", df_name))
  }

  if (any(grepl(" ", parameters$Variable, fixed = TRUE))) {
    chk::err(sprintf("spaces are not allowed in the names present in the `Variable` column in %s", df_name))
  }

  if (n_proxies > 0 && any(paste0("p", seq_len(n_proxies)) %in% parameters$Variable)) {
    chk::err(sprintf("to avoid name clashes with the proxy variable%%s, please avoid using %s in the `Variable` column in %s",
                     word_list(paste0("p", seq_len(n_proxies)), and.or = "or", quotes = 2),
                     df_name),
             n = n_proxies)
  }

  #Type
  parameters$Type <- as.character(parameters$Type)
  if (anyNA(parameters$Type) || any(parameters$Type == "")) {
    chk::err(sprintf("missing and empty values are not allowed in the `Type` column in %s", df_name))
  }

  if (!all(parameters$Type %in% c("binary", "continuous", "count"))) {
    chk::err(sprintf("all values in the `Type` column in %s must be \"binary\", \"continuous\", or \"count\"",
                     df_name))
  }

  binary <- parameters$Type == "binary"
  continuous <- parameters$Type == "continuous"
  count <- parameters$Type == "count"

  #prevalence
  parameters$prevalence <- suppressWarnings(as.numeric(parameters$prevalence))
  if (any(binary)) {
    if (anyNA(parameters$prevalence[binary])) {
      chk::err(sprintf("values in the `prevalence` column of %s must be supplied for all binary variables present",
                       df_name))
    }

    if (any(parameters$prevalence[binary] <= 0) || any(parameters$prevalence[binary] >= 1)) {
      chk::err(sprintf("all values in the `prevalence` column of %s must be between 0 and 1",
                       df_name))
    }
  }

  if (any(continuous)) {
    if (!all(is.na(parameters$prevalence[continuous]))) {
      chk::wrn(sprintf("values in the `prevalence` column of %s for continuous variables will be ignored",
                       df_name))
    }

    is.na(parameters$prevalence[continuous]) <- TRUE
  }

  if (any(count)) {
    if (!all(is.na(parameters$prevalence[count]))) {
      chk::wrn(sprintf("values in the `prevalence` column of %s for count variables will be ignored",
                       df_name))
    }

    is.na(parameters$prevalence[count]) <- TRUE
  }

  #mean
  parameters$mean <- suppressWarnings(as.numeric(parameters$mean))
  if (any(continuous)) {
    if (anyNA(parameters$mean[continuous])) {
      chk::err(sprintf("values in the `mean` column of %s must be supplied for all continuous variables present",
                       df_name))
    }
  }

  if (any(count)) {
    if (anyNA(parameters$mean[count])) {
      chk::err(sprintf("values in the `mean` column of %s must be supplied for all count variables present",
                       df_name))
    }

    if (any(parameters$mean[count] <= 0)) {
      chk::err(sprintf("all values in the `mean` column of %s must be greater than 0 for count variables",
                       df_name))
    }
  }

  if (any(binary)) {
    if (!all(is.na(parameters$mean[binary]))) {
      chk::wrn(sprintf("values in the mean column of %s for binary variables will be ignored",
                       df_name))
    }

    is.na(parameters$mean[binary]) <- TRUE
  }

  #sd
  parameters$sd <- suppressWarnings(as.numeric(parameters$sd))
  if (any(continuous)) {
    if (anyNA(parameters$sd[continuous])) {
      chk::err(sprintf("values in the `sd` column of %s must be supplied for all continuous variables present",
                       df_name))
    }

    if (any(parameters$sd[continuous] <= 0)) {
      chk::err(sprintf("all values in the `sd` column of %s must be greater than 0",
                       df_name))
    }
  }

  if (any(count)) {
    if (!all(is.na(parameters$sd[count]))) {
      chk::wrn(sprintf("values in the `sd` column of %s for count variables will be ignored",
                       df_name))
    }

    is.na(parameters$sd[count]) <- TRUE
  }

  if (any(binary)) {
    if (!all(is.na(parameters$sd[binary]))) {
      chk::wrn(sprintf("values in the `sd` column of %s for binary variables will be ignored",
                       df_name))
    }

    is.na(parameters$sd[binary]) <- TRUE
  }

  #coeff_treatment_model/coeff_outcome_model
  for (i in c("coeff_treatment_model", "coeff_outcome_model")) {
    nas <- is.na(parameters[[i]])
    parameters[[i]] <- suppressWarnings(as.numeric(parameters[[i]]))

    if (anyNA(parameters[[i]][!nas])) {
      chk::err(sprintf("all values in the `%s` column of %s must be numeric or empty",
                       i, df_name))
    }

    parameters[[i]][nas] <- 0

    if (parameters[[i]][parameters$Variable == unmeasured_conf] == 0) {
      chk::err(sprintf("the value of `%s` in %s cannot be empty or zero for the variable named in `unmeasured_conf` (\"%s\")",
                       i, df_name, unmeasured_conf))
    }
  }

  if (any(parameters$coeff_treatment_model == 0 & parameters$coeff_outcome_model == 0)) {
    chk::err("every variable in %s must have a nonzero value of at least one of `coeff_treatment_model` and `coeff_outcome_model`",
             df_name)
  }

  parameters
}
