simBA_de_novo <- function(iterations, parameter_file_path, size, treatment_prevalence, treatment_coeff,
                          outcome_prevalence, dist='E',
                          unmeasured_conf, n_proxies, proxy_type, corr){

  ### function to iterate the simulation runs

  simBA_runs <- function (iterations, parameter_file_path, size, treatment_prevalence,
                          treatment_coeff, outcome_prevalence, dist,
                          unmeasured_conf, n_proxies, proxy_type, corr){
    HRs = data.frame()
    bal = data.frame()

    for(i in 1:iterations)
    {
      try(res_i <- simulate_analyze_u(parameter_file_path, size, treatment_prevalence, treatment_coeff, outcome_prevalence, dist,
                                      unmeasured_conf, n_proxies, proxy_type, corr))

      HRs <- rbind(HRs, res_i[[2]])
      bal <- rbind(bal, res_i[[3]])
    }

    return(list(bal, HRs))

  }


  ### run the function and save the results

  runs_results <- capture.output({run_results <- simBA_runs(iterations, parameter_file_path, size, treatment_prevalence, treatment_coeff, outcome_prevalence, dist,
                                                            unmeasured_conf, n_proxies, proxy_type, corr)})
  ## prepare the the SMD table to output

  balance <- run_results[[1]]

  balance_table <- balance[(balance$SMD != ''),]

  balance_table$SMD <- as.numeric(balance_table$SMD)

  balance_table$SMD[is.na(balance_table$SMD)] <- 0

  SMD_table <- balance_table %>% group_by(variable, adjustment) %>%
    summarise(across(c(SMD),mean),
              .groups = 'drop') %>%
    as.data.frame()

  ## prepare the HR table to output

  results_n <- run_results[[2]]

  HR_table <- results_n %>% group_by(adjustment) %>%
    summarise(across(c(HR),mean),
              .groups = 'drop') %>%
    as.data.frame()

  HR_table$HR <- exp(HR_table$HR)

  return(list(SMD_table, HR_table))

}

