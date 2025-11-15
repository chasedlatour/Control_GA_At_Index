###########################################################################
#### FILE: data_generation.R
#### PURPOSE: The functions created in this files are to be used to create 
#### a simulated dataset of pregnancies.
###########################################################################



######################################################################
# FUNCTION: generate_pregs
# PURPOSE: The purpose of this function is to generate pregnancies for
# the simulation study 
# INPUTS: 
# - sim_id = identifier of the replicate number of the cohort 
# - n = Number of pregnancies to generate
# - z_prob = baseline probability of Z
#
# 
sim_id <- 1
n <- 200
prob_z <- 0.5
p_ptb_z0<- 0.1 # probability of delivery by 37 weeks gestation (i.e., preterm birth)
p_ptb_z1<- 0.2
first_ga <- 15
weeks <- seq(15, 22, by=1)
prob_enter <- seq(0.125, 8)





generate_pregs <- function(sim_id, n, prob_z, p_ptb_z0, p_ptb_z1,
                           prob_enter, first_ga, weeks,
                           shape_cens_z0, scale_cens_z0,
                           shape_cens_z1, scale_cens_z1){
  
  # Record key identifiers
  n_sim <- sim_id   # identifier associated with the simulated dataset
  id <- 1:n         # person-level identifier
  
  # Determine value of binary Z
  z <- rbinom(n=n, size=1, prob=prob_z)
  
  
  ### GENERATE OUTCOMES STARTING AT 16 WEEKS OF GESTATION
  
  ## Z == 0
  # Create parameters for preterm birth generation
  p_event <- p_ptb_z0 # probability of an event (preterm birth) 
  birthtau <- (37-first_ga)/52 #Maximum pregnancy duration in years
  alpha <- 0.4 # weibull shape for time from 16 weeks of gestation to outcome
  lambda <- (-log(1-p_event))^(1/alpha)/birthtau   # weibull scale for time to PTB
  h <- exp(log(lambda)) 
  # Generate the time to delivery
  time_to_event_z0 <- rweibull(n, shape=alpha, scale=1/h)        # time to delivery
  
  # Create an indicator for PTB in a setting without censoring
  ptb_nocens_z0 <- ifelse(time_to_event_z0 > birthtau, 0, 1)
  tte_ptb_nocens_z0 <- ifelse(ptb_nocens_z0 == 1, time_to_event_z0, birthtau)
  
  ## Z == 1
  
  # Create parameters for preterm birth generation
  p_event <- p_ptb_z1 # probability of an event (preterm birth) 
  birthtau <- (37-first_ga)/52 #Maximum pregnancy duration in years
  alpha <- 0.4 # weibull shape for time from 16 weeks of gestation to outcome
  lambda <- (-log(1-p_event))^(1/alpha)/birthtau   # weibull scale for time to PTB
  h <- exp(log(lambda)) 
  # Generate the time to delivery
  time_to_event_z1 <- rweibull(n, shape=alpha, scale=1/h)        # time to delivery
  
  # Create an indicator for PTB in a setting without censoring
  ptb_nocens_z1 <- ifelse(time_to_event_z1 > birthtau, 0, 1)
  tte_ptb_nocens_z1 <- ifelse(ptb_nocens_z1 == 1, time_to_event_z1, birthtau)
  
  ## Indicators based on observed value of Z
  ptb_nocens <- ifelse(z==1, ptb_nocens_z1, ptb_nocens_z0)
  tte_ptb_nocens <- ifelse(z==1, tte_ptb_nocens_z1, tte_ptb_nocens_z0)
  
 
  ### GENERATE THEIR ACTUAL GESTATIONAL WEEK AT ENROLLMENT
  
  mat <- rmultinom(n=n, size=1, prob_enter)
  ga_entry <- weeks[max.col(t(mat))]
  
  ### GENERATE LTFU AT THE START OF ENROLLMENT FOR EACH PERSON
  
  # Generate the time to censoring under Z=0
  time_to_censor_z0 <- rweibull(n, shape=shape_cens_z0, scale=scale_cens_z0)
  # Generate the time to censoring under Z=1
  time_to_censor_z1 <- rweibull(n, shape=shape_cens_z1, scale=scale_cens_z1)
  
  # Observed time-to-censor
  time_to_censor <- ifelse(z==1, time_to_censor_z1, time_to_censor_z0)
  
  ### NOW DETERMINE THEIR OBSERVED OUTCOME, APPLYING TIME TO CENSORING AFTER ENROLLMENT
  time_to_censor_ga <- (ga_entry - first_ga)/52 + time_to_censor
  # Observed event times and types if all enrolled at 16 weeks gestation
  tte <- pmin(tte_ptb_nocens, time_to_censor_ga, birthtau)
  t_weeks <- first_ga + (tte* 52)   # Gestational age at end of follow-up
  outcome <- ifelse(ptb_nocens == 1 & tte == tte_ptb_nocens,
                    1,
                    ifelse(
                      ptb_nocens == 0 & tte == birthtau,
                      2,
                      0
                    ))
  
 
  #### Make the base dataset

  data <- tibble(n_sim, id, ga_entry, z, ptb_nocens, outcome, t_weeks) %>% 
    # Calculate time to event from ga_entry
    mutate(tte = t_weeks - ga_entry) %>% 
    # Exclude those individuals who experienced a delivery before trial entry
    filter(tte >= 0)
  
  return(data)
}


test <- generate_pregs(sim_id = 1, n=20000, prob_z = 0.5, 
                       p_ptb_z0 = 0.15, p_ptb_z1 = 0.3, prob_enter = rep(0.125, 8),
                       first_ga = 15, weeks = seq(15,22, by=1),
                       shape_cens_z0 = 3, scale_cens_z0 = 0.5,
                       shape_cens_z1 = 5, scale_cens_z1 = 0.6)


prop.table(table(test$z, test$outcome))
prop.table(table(test$ga_entry, test$outcome))
censor <- subset(test, outcome == 0)
prop.table(table(test$ga_entry, test$z))


# Estimate the probability of preterm birth without censoring

mean(test$ptb_nocens)

# Estimate the probability of preterm birth using an AJ estimator, not adjusting for Z
mod <- summary(survfit(Surv(tte, factor(outcome)) ~ z1, data=test))


# Estimate the probability of preterm birth using an AJ estimator, stratify on Z and then
# standardize to the distribution of Z at baseline. Censoring random within strata of Z

mod <- summary(survfit(Surv(tte, factor(outcome)) ~ z, data=test))
# Now get the cumulative incidence overall
prob <- data.frame(z_strat = mod$strata, r = mod$pstate[,2]) %>% 
  mutate(
    z = case_when(z_strat == "z=0" ~ 0,
                  z_strat == "z=1" ~ 1)
  ) %>% 
  group_by(z) %>% 
  slice_tail(n=1) 

# Derive the probabilities of each value of Z at baseline
props <- test %>% 
  group_by(z) %>% 
  summarize(prop=n()/nrow(test))

# Merge them together
merge <- left_join(prob, props, by = c("z"="z"))

standardized <- merge %>% 
  mutate(std_risk = r*prop) %>% 
  ungroup() %>% 
  summarize(risk = sum(std_risk))
standardized




# Now estimate risks within each GA stratum and standardize

mod <- summary(survfit(Surv(tte, factor(outcome)) ~ ga_entry, data=test))
strat_risk <- tibble(ga = mod$strata,r = mod$pstate[,2]) %>% 
  mutate(
    ga_num = case_when(ga == "ga_entry=15" ~ 15,
                       ga == "ga_entry=16" ~ 16,
                       ga == "ga_entry=17" ~ 17,
                       ga == "ga_entry=18" ~ 18,
                       ga == "ga_entry=19" ~ 19,
                       ga == "ga_entry=20" ~ 20,
                       ga == "ga_entry=21" ~ 21,
                       ga == "ga_entry=22" ~ 22)
  ) %>% 
  group_by(ga_num) %>% 
  slice_tail(n=1) 

# Derive the probabilities of each gestational week at baseline
props <- test %>% 
  group_by(ga_entry) %>% 
  summarize(prop = n()/nrow(test))

# Merge them together
merge <- left_join(strat_risk, props, by = c("ga_num"="ga_entry"))

standardized <- merge %>% 
  mutate(std_risk = r*prop) %>% 
  ungroup() %>% 
  summarize(risk = sum(std_risk))
  

# output side by side

mean(test$ptb_nocens)
prob
standardized






















#OLD


# WITH CNFOUNDING AND TREATMENT
# sim_id <- 1 
# n <- 200
# z_prob <- 0.5
# prob_trt_z01 <- c(0.4, 0.6)
# prob_enter_z0 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
# prob_enter_z1 <- c(0.3, 0.3, 0.2, 0.1, 0.1)
# 
# generate_pregs <- function(sim_id, n, z_prob, prob_trt_z01,
#                            prob_enter_z0, prob_enter_z1){
#   
#   p_entry <- rbind(prob_enter_z0, prob_enter_z1)
#   
#   # Create the baseline dataset
#   data <- tibble(
#     n_sim = sim_id,
#     id = 1:n,
#     
#     # Determine each pregnancy's Z value
#     z = rbinom(n=n, size=1, prob = z_prob),
#     
#     # Determine the probability of treatment dependent upon the observed value of Z
#     p_trt = prob_trt_z01[z+1],
#     # Now assign treatment
#     trt = rbinom(n=n, size=1, prob=p_trt),
#     
#     # Record probabilities for study entry
#     p_entry = ifelse(z == 1,
#                      list(prob_enter_z1),
#                      list(prob_enter_z0)
#     ),
#     
#     # Choose the week of trial entry
#     wk_entry = map_int(p_entry, ~ {
#       draw <- rmultinom(1, size = 1, prob = .x)
#       which(draw == 1)
#     }),
#     ga_entry = case_when(wk_entry == 1 ~ 16,
#                          wk_entry == 2 ~ 17,
#                          wk_entry == 3 ~ 18,
#                          wk_entry == 4 ~ 19,
#                          wk_entry == 5 ~ 20)
#   ) 
# }

