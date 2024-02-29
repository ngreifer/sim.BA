library(readxl)
library(msm)
library(simstudy)
library(survival)
library(MatchIt)
library(reshape)
library(tableone)
library(ggplot2)
library(ggthemes)
library(dplyr)

### Function 1- simulate the data #####

simulate_analyze_u <- function (parameter_file_path, size, treatment_prevalence, treatment_coeff, outcome_prevalence, dist='E',
                                unmeasured_conf, n_proxies, proxy_type, corr){

  # Read the Excel file
  df <- read_excel(parameter_file_path)


  # Task 1: generate covariates with pre-specified prevalence if binary and mean/sd if continuous #####

  # split out binary and continuous input variables

  df1 <- df[df$Type == "binary", ]
  df1$obs <-  1:nrow(df1) #useful for do loops to generate these variables
  df2 <- df[df$Type != "binary", ]
  df2$obs <-  1:nrow(df2) #useful for do loops to generate these variables

  generate_binary <- function(prev){
    values <- data.frame(rbinom(size, 1, prev))
    return(values)
  }

  out_binary <- c()
  for (i in 1:max(df1$obs)) {
    (out_binary[i] <- generate_binary(prev = as.numeric(df1[i, "prevalence"])))
  }

  binaries <- data.frame(out_binary)

  #attach names to the generated variables based on user specified 'variable' column

  names(binaries) <- paste0(df1$Variable)


  generate_cont <- function(mean, sd){
    values <- data.frame(rnorm(size, mean, sd))
    return(values)
  }

  out_cont <- c()
  for (i in 1:max(df2$obs)) {
    (out_cont[i] <- generate_cont(mean = as.numeric(df2[i, "mean"]), sd= as.numeric(df2[i, "sd"])))
  }

  continuous <- data.frame(out_cont)

  #attach names to the generated variables based on user specified 'variable' column

  names(continuous) <- paste0(df2$Variable)

  simdata <- cbind(binaries, continuous)


  # Task 2: write out treatment and outcome models based on user specified coefficients #####

  ## binary treatment model, extract supplied coefficients to throw in the equation

  df_treatment_model_subset <- subset(df, !is.na(coeff_treatment_model))
  df_treatment_model_subset$v1 <- paste("simdata", df_treatment_model_subset$Variable, sep="$")

  df_treatment_model_subset$treatment_model_terms <- paste(df_treatment_model_subset$v1, df_treatment_model_subset$coeff_treatment_model, sep = "*")
  exp_term=paste(df_treatment_model_subset$treatment_model_terms, collapse = '+')
  simdata$exp_term <- eval(parse(text=exp_term))

  # a loop to match user desired prevalence for the treatment

  intercept_init <- -0 # initiate with intercept at 0, this will be updated in the loop below to match the desired prevalence

  simdata$trueps_init <- (1 + exp( -(intercept_init + simdata$exp_term)))^-1 #initiate true PS calc with 0 intercept

  # probability of exposure random number betw 0 and 1, estimate the generated treatment prevalence

  prob.exposure <- runif(size)
  simdata$treatment_init <- ifelse(simdata$trueps_init > prob.exposure, 1, 0)
  generated_treatment_prev <- sum(simdata$treatment_init)/size

  #take a 1000 tries varying intercept by 0.02 each time in either direction depending on the inequality between desired and generated prevalence

  for (i in 1:1000) {
    if (treatment_prevalence <= generated_treatment_prev) { #lower the intercept if desired prev is less than or equal to generated
      intercept <- intercept_init - 0.02 * i
      simdata$trueps <- (1 + exp(-(intercept + simdata$exp_term))) ^ -1
      simdata$treatment <- ifelse(simdata$trueps > prob.exposure, 1, 0)
      generated_treatment_prev <- sum(simdata$treatment) / size
    } else { #increase the intercept if desired prev is less than or equal to generated
      intercept <- intercept_init + 0.02 * i
      simdata$trueps <- (1 + exp(-(intercept + simdata$exp_term))) ^ -1
      simdata$treatment <- ifelse(simdata$trueps > prob.exposure, 1, 0)
      generated_treatment_prev <- sum(simdata$treatment) / size
    }
    if ((generated_treatment_prev-0.002) <= (treatment_prevalence) & (treatment_prevalence) <= (generated_treatment_prev+0.002) ) {   ## stop when generated and desired prevalences match
      break
    }
    #print(i)
    #print(generated_treatment_prev)
  }

  #save the simulated treatment prevalence

  sim_exp_prevalence <- sum(simdata$treatment)/size

  #sim_exp_prevalence

  #exclude the helpers

  simdata <- simdata[, !(names(simdata) %in% c("trueps_init", "treatment_init", "exp_term"))]


  ## time to event outcome, extract supplied coefficients to throw in the equation

  df_outcome_model_subset <- subset(df, !is.na(coeff_outcome_model))
  df_outcome_model_subset$v1 <- paste("simdata", df_outcome_model_subset$Variable, sep="$")
  df_outcome_model_subset$outcome_model_terms <- paste(df_outcome_model_subset$v1, df_outcome_model_subset$coeff_outcome_model, sep = "*")
  linear_pred =paste(df_outcome_model_subset$outcome_model_terms, collapse = '+')
  simdata$linear_pred <- eval(parse(text=linear_pred))

  #need to update linear predictor with user specified treatment effect
  simdata$linear_pred_w_treatment <- simdata$linear_pred + simdata$treatment*treatment_coeff


  #building the outcome time-to-event model

  Uevent <- runif (size,0,1) #Random variable distributed uniformly, meant to represent survival function

  simdata$covariate.vector <- exp (simdata$linear_pred_w_treatment)

  scale.event <- 0.0001 # scale parameter for exponentially distributed event time

  #now we draw an event time and a censoring time, the minimum  of these is "observed" and we record whether it was the event time or the censoring time

  if (dist=='E') {
    simdata$true_time <- (-log(Uevent))/(scale.event*simdata$covariate.vector)
  } else if (dist=='W') {
    simdata$true_time <- ((-log(Uevent))/(scale.event*simdata$covariate.vector))^0.5
  }

  #generate administrative censoring times to generate outcome data

  simdata$cens_time <- quantile(simdata$true_time, c(outcome_prevalence)) #censor at nth percentile of the generated event time to generate desired outcome prevalence

  simdata$survt = pmin(simdata$cens_time, simdata$true_time)  #observed time is min of censored and true

  simdata$event = simdata$survt==simdata$true_time   # set to 1 if event is observed

  table(simdata$event)

  #exclude the helpers

  simdata <- simdata[, !(names(simdata) %in% c("true_time", "cens_time", "linear_pred", "linear_pred_w_treatment", "covariate.vector"))]


  ## TASK 3- insert proxies ########

  if (n_proxies>0){

    out_proxies <- c()


    if (proxy_type== 'binary'){

      add_binary_proxy = function(x, rho) {

        #Calculate the mean of x
        x_bar = mean(x)

        #Calculate the expected value of xy under the assumption that x and y are correlated with correlation coefficient rho
        xy_bar = rho * x_bar + (1 - rho) * x_bar ^ 2

        #Calculate the number of 1s in x that need to be flipped to achieve the desired correlation coefficient
        toflip = sum(x == 1) - round(length(x) * xy_bar)

        #Create a copy of x
        y = x

        #Randomly flip the values of x at the specified indexes to achieve the desired correlation coefficient
        y[sample(which(x == 0), toflip)] = 1
        y[sample(which(x == 1), toflip)] = 0

        proxy <- data.frame(y)
        #Return the modified vector
        return(proxy)
      }

      for (i in 1:n_proxies) {
        (out_proxies[i] <- add_binary_proxy(simdata[,(unmeasured_conf)], corr))
      }
    }else if (proxy_type== 'continuous') {

      add_continuous_proxy <- function(x, rho)
      {
        y <- (rho * (x - mean(x)))/sqrt(var(x)) + sqrt(1 - rho^2) * rnorm(length(x))

        return(y)
      }

      for (i in 1:n_proxies) {
        (out_proxies[i] <- add_continuous_proxy(simdata[,(unmeasured_conf)], corr))
      }
    }

    proxies <- data.frame(out_proxies)

    for (i in 1:n_proxies) {
      # Get the name of the i'th proxy
      name = paste("p", i, sep = "")

      # Assign the name to the i'th proxy
      names(proxies)[[i]] = name
    }


    simdata_p <- cbind(simdata, proxies)

  }else if(n_proxies==0){
    simdata_p <- simdata
  }

  ## TASK 4- analysis ###########

  ## model 1- crude

  crude_model <- coxph(formula = Surv(simdata$survt, simdata$event) ~ simdata$treatment, method = "breslow")
  HR= summary(crude_model)$coefficients[1,1]
  SE= summary(crude_model)$coefficients[1,3]
  crude_results <- data.frame(HR, SE)
  crude_results$adjustment <- 'crude'

  ## model 2- adjusted for all measured confounding (L1)

  #exclude the ones we dont want to adjust for

  simdata_adjustment <- simdata[,-which(names(simdata) %in% c("trueps", "treatment", "survt", "event", unmeasured_conf))]

  # write the formula

  adjust_variables <- colnames(simdata_adjustment[c(1:ncol(simdata_adjustment))])
  f <- as.formula(
    paste("treatment",
          paste(adjust_variables, collapse = " + "),
          sep = " ~ "))

  m.out <- matchit(f, data=simdata_p,
                   method="nearest",
                   caliper=0.2,
                   replace=FALSE)

  matched_L1 <- match.data(m.out) #save PS matched data

  L1_model <- coxph(formula = Surv(matched_L1$survt, matched_L1$event) ~ matched_L1$treatment, method = "breslow")
  HR= summary(L1_model)$coefficients[1,1]
  SE= summary(L1_model)$coefficients[1,3]
  L1_results <- data.frame(HR, SE)
  L1_results$adjustment <- 'L1'

  ## model 3- adjusted for all measured confounding (L1) plus the proxies

  if (n_proxies>0){
    simdata_adjustment_proxy <- cbind(simdata_adjustment, proxies)


    # write the formula

    adjust_variables_proxy <- colnames(simdata_adjustment_proxy[c(1:ncol(simdata_adjustment_proxy))])
    f1 <- as.formula(
      paste("treatment",
            paste(adjust_variables_proxy, collapse = " + "),
            sep = " ~ "))

    m.out1 <- matchit(f1, data=simdata_p,
                      method="nearest",
                      caliper=0.2,
                      replace=FALSE)

    matched_L2 <- match.data(m.out1) #save PS matched data

    L2_model <- coxph(formula = Surv(matched_L2$survt, matched_L2$event) ~ matched_L2$treatment, method = "breslow")
    HR= summary(L2_model)$coefficients[1,1]
    SE= summary(L2_model)$coefficients[1,3]
    L2_results <- data.frame(HR, SE)
    L2_results$adjustment <- 'L2'

    all_HRs <- rbind(crude_results, L1_results, L2_results)

  } else if (n_proxies==0){
    all_HRs <- rbind(crude_results, L1_results)
  }



  #balance measures
  if (n_proxies>0){
    balance_vars <- paste(names(proxies))

    balance_crude <- CreateTableOne(data = simdata_p, strata = "treatment", vars = c(unmeasured_conf, balance_vars))
    table_balance_crude <- as.data.frame(print(balance_crude, smd = T))
    table_balance_crude$adjustment <- 'crude'
    table_balance_crude$variable <- row.names(table_balance_crude)
    rownames(table_balance_crude)<-NULL

    balance_L1<- CreateTableOne(data = matched_L1, strata = "treatment", vars = c(unmeasured_conf, balance_vars))
    table_balance_L1 <- as.data.frame(print(balance_L1, smd = T))
    table_balance_L1$adjustment <- 'L1'
    table_balance_L1$variable <- row.names(table_balance_L1)
    rownames(table_balance_L1)<-NULL

    balance_L2<- CreateTableOne(data = matched_L2, strata = "treatment", vars = c(unmeasured_conf, balance_vars))
    table_balance_L2 <- as.data.frame(print(balance_L2, smd = T))
    table_balance_L2$adjustment <- 'L2'
    table_balance_L2$variable <- row.names(table_balance_L2)
    rownames(table_balance_L2)<-NULL

    all_balance_tables <- rbind(table_balance_crude, table_balance_L1, table_balance_L2)

  } else if (n_proxies==0){
    balance_crude <- CreateTableOne(data = simdata_p, strata = "treatment", vars = c(unmeasured_conf))
    table_balance_crude <- as.data.frame(print(balance_crude, smd = T))
    table_balance_crude$adjustment <- 'crude'
    table_balance_crude$variable <- row.names(table_balance_crude)
    rownames(table_balance_crude)<-NULL

    balance_L1<- CreateTableOne(data = matched_L1, strata = "treatment", vars = c(unmeasured_conf))
    table_balance_L1 <- as.data.frame(print(balance_L1, smd = T))
    table_balance_L1$adjustment <- 'L1'
    table_balance_L1$variable <- row.names(table_balance_L1)
    rownames(table_balance_L1)<-NULL

    all_balance_tables <- rbind(table_balance_crude, table_balance_L1)
  }


  if (n_proxies>0){
    return(list(simdata_p, all_HRs, all_balance_tables))
  }else if (n_proxies==0){
    return(list(simdata, all_HRs, all_balance_tables))
  }


}



