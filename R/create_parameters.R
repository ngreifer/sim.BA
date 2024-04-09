#' Create parameters for data-generating models
#'
#' `create_parameters()` facilitates creation of the `parameters` input to [simBA()]. This input contains information required to generate simulated dataset that will be analyzed within the simulation.
#'
#' @param nbinary,ncontinuous,ncount the number of binary, continuous, and count confounders to include, respectively. Default is 0.
#' @param unmeasured_conf the name of the unmeasured confounder. Default is `"u1"`.
#' @param unmeasured_type the type of variable for the unmeasured confounder. Allowable options include `"binary"`, `"continuous"`, and `"count"`. Default is `"binary"`; abbreviations allowed.
#' @param file optional; a string containing a path to a .csv or .xslx file where the output will be written to. If `NULL` (the default), no file will be written.
#'
#' @returns
#' A data.frame containing a skeleton of the parameter values, which must be filled in manually by the user. See *The `parameters` input object* section in the [simBA()] documentation for which columns will be present in the output. An additional column, Description, will also be produced, but it is not necessary to fill it in.
#'
#' @examples
#' parameters <- create_parameters(nbinary = 6,
#'                                 ncontinuous = 2,
#'                                 ncount = 1,
#'                                 unmeasured_conf = "u1",
#'                                 unmeasured_type = "continuous")
#'
#' parameters


#' @export
create_parameters <- function(nbinary = 0, ncontinuous = 0, ncount = 0, unmeasured_conf = "u1",
                              unmeasured_type = "binary", file = NULL) {

  chk::chk_count(nbinary)
  chk::chk_count(ncontinuous)
  chk::chk_count(ncount)

  chk::chk_string(unmeasured_conf)
  chk::chk_string(unmeasured_type)
  unmeasured_type <- tolower(unmeasured_type)
  unmeasured_type <- match_arg(unmeasured_type, c("binary", "continuous", "count"))

  parameter_names <- c("Variable", "Description", "Type", "prevalence", "mean", "sd",
                      "coeff_treatment_model", "coeff_outcome_model")

  nvars <- nbinary + ncontinuous + ncount

  var_names <- paste0("c", seq_len(nvars), recycle0 = TRUE)

  if (unmeasured_conf %in% var_names) {
    chk::err("please provide another name to `unmeasured_conf`")
  }

  parameters <- as.data.frame(matrix(nrow = nvars + 1, ncol = length(parameter_names),
                                     dimnames = list(NULL, parameter_names)))

  parameters$Variable <- c(var_names, unmeasured_conf)

  parameters$Description <- ""
  parameters$Description[parameters$Variable == unmeasured_conf] <- "Unmeasured"

  parameters$Type <- c(rep(c("binary", "continuous", "count"), c(nbinary, ncontinuous, ncount)),
                       unmeasured_type)

  if (!is.null(file)) {
    chk::chk_string(file)
    if (endsWith(file, ".csv")) {
      utils::write.csv(parameters, file, row.names = FALSE)
    }
    else if (endsWith(file, ".xlsx")) {
      rlang::check_installed("openxlsx")

      openxlsx::write.xlsx(parameters, file)
    }
    else {
      chk::err("if supplied, `file` must refer to a .csv or .xlsx file")
    }

    return(invisible(parameters))
  }

  parameters
}