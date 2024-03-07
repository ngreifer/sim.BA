#' Plot the results of a simulation
#'
#' `plot()` plots the output of a call to `simBA()`. The plot can contain either the estimated hazard ratios or standardized mean differences across simulations, each in a set of box plots.
#'
#' @param x a `simBA` object; the output of a call to [simBA()].
#' @param type the type of plot to produce; allowable options include `"balance"` (default) and `"hr"`. Abbreviations allowed.
#' @param ... further arguments passed to [ggplot2::geom_boxplot()].
#'
#' @returns
#' A `ggplot` object, which can be modified using `ggplot2` syntax.
#'
#' @details
#' The balance plot plots absolute standardized mean differences. Vertical lines are placed at 0 (solid) and .1 (dashed). The hazard ratio (HR) plot plots hazard ratios on a log scale for the x-axis. Vertical lines are placed at 1 (solid) and the true marginal HR (dashed).
#'
#' @seealso [simBA()] for performing the simulation.
#'
#' @examples
#' # See help("simBA") for examples.
#'

#' @exportS3Method plot simBA
plot.simBA <- function(x, type = "balance", ...) {

  chk::chk_string(type)
  type <- tolower(type)
  type <- match_arg(type, c("balance", "hr"))
  adj <- attr(x, "adj")

  if (type == "balance") {
    bal <- do.call("rbind", lapply(x$sim_out, `[[`, "all_balance_tables"))

    proxy_vars <- sort(setdiff(unique(bal$variable), attr(x, "unmeasured_conf")))
    new_proxy_vars <- sub("p", "Proxy ", proxy_vars, fixed = TRUE)

    ## prepare the the SMD table to output
    bal$SMD[is.na(bal$SMD)] <- 0

    bal$variable <- factor(bal$variable,
                            levels = c(attr(x, "unmeasured_conf"), rev(proxy_vars)),
                            labels = c("Unmeasured\nConfounder", rev(new_proxy_vars)))
    bal$adjustment <- factor(bal$adjustment, levels = c("crude", "L1", "L2"),
                             labels = c("Unadjusted", paste("Level 1", firstup(adj)), paste("Level 2", firstup(adj))))

    bal$scenario <- factor(paste(bal$variable, bal$adjustment),
                           nmax = nlevels(bal$variable) * nlevels(bal$adjustment ))

    p <- ggplot(bal) +
      geom_boxplot(aes(x = abs(.data$SMD),
                       y = .data$variable,
                       fill = .data$variable),
                   ...) +
      geom_vline(xintercept = 0) +
      geom_vline(xintercept = .1, linetype = "dashed") +
      geom_vline(xintercept = mean(abs(bal$SMD[bal$variable == "Unmeasured\nConfounder" &
                                                   bal$adjustment == "Unadjusted"])),
                 linetype = "dotted") +
      facet_grid(~adjustment, switch = "y") +
      theme_bw() +
      theme(strip.placement = "outside") +
      guides(fill = "none") +
      labs(x = "Absolute Standardized\nMean Difference", y = NULL)

    p

  }
  else { #type == "hr"
    est <- do.call("rbind", lapply(x$sim_out, `[[`, "all_HRs"))

    est <- est[est$adjustment != "true",]

    est$adjustment <- factor(est$adjustment, levels = c("L2", "L1", "crude"),
                             labels = c(paste("Level 2", firstup(adj)), paste("Level 1", firstup(adj)), "Crude"))

    p <- ggplot(est) +
      geom_boxplot(aes(x = exp(.data$HR),
                       y = .data$adjustment,
                       fill = .data$adjustment),
                   ...) +
      geom_vline(xintercept = 1) +
      geom_vline(xintercept = x$HR_table$HR[x$HR_table$adjustment == "true"],
                 linetype = "dashed") +
      scale_x_log10() +
      theme_bw() +
      guides(fill = "none") +
      labs(x = "Hazard Ratio", y = NULL)

    p
  }
}
