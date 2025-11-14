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

sim_id <- 1 
n <- 200
prob_enter <- c(0.2, 0.2, 0.2, 0.2, 0.2)
p_ptb<- 0.2 # probability of delivery by 37 weeks gestation (i.e., preterm birth)
study_tau <- 21/52# last follow up time for preterm birth (in years)
cenrate <- 0.1 # How much censoring



generate_pregs <- function(sim_id, n, prob_enter, p_ptb){
  
  # Record key identifiers
  n_sim <- sim_id   # identifier associated with the simulated dataset
  id <- 1:n         # person-level identifier
  
  
  ### GENERATE OUTCOMES UNDER THE ASSUMPTION THAT ALL ENROLL AT 16 WEEKS GESTATION
  
  # Create parameters for preterm birth generation
  p_event <- p_ptb # probability of an event (preterm birth) 
  birthtau <- 21/52 #Maximum pregnancy duration in years
  alpha <- 0.5 # weibull shape for time from 16 weeks of gestation to outcome
  lambda <- (-log(1-p_event))^(1/alpha)/birthtau   # weibull scale for time to PTB
  h <- exp(log(lambda)) 
  
  # Generate the time to delivery
  time_to_event <- rweibull(n, shape=alpha, scale=1/h)        # time to delivery
  
  # Create an indicator for PTB in a setting without censoring
  ptb_nocens <- ifelse(time_to_event > birthtau, 0, 1)
  tte_ptb_nocens <- ifelse(ptb_nocens == 1, time_to_event, birthtau)
  
  # Create parameters for censoring generation
  p_censor <- cenrate
  alpha_cens <- 0.22 # shape parameter for time to censoring
  lambda_cens <- (-log(1-p_censor))^(1/alpha_cens)/birthtau # weibull scale for time to censoring
  h_cens <- exp(log(lambda_cens))
  
  # Generate the time to censoring
  time_to_censor <- rweibull(n, shape=alpha_cens, scale=1/h_cens)
  
  # # Create time to censoring
  # u <- runif(n, 0, 1)
  # c <- runif(n=n, min=0, max=birthtau-0.001)*(u<cenrate) + birthtau*(u>=cenrate)  

  # Observed event times and types if all enrolled at 16 weeks gestation
  tte <- pmin(tte_ptb_nocens, time_to_censor, birthtau)
  t_weeks <- 16+(tte* 52)   # Gestational age at end of follow-up
  outcome <- ifelse(ptb_nocens == 1 & tte == tte_ptb_nocens,
                    1,
                    ifelse(
                      ptb_nocens == 0 & tte == birthtau,
                      2,
                      0
                    ))
  
  # Now convert everything to weeks, under the assumption that there are 52 weeks in a year
  
  
  ### GENERATE THEIR ACTUAL GESTATIONAL WEEK AT ENROLLMENT
  
  mat <- rmultinom(n=n, size=1, c(0.2, 0.2, 0.2, 0.2, 0.2))
  wk_entry <- max.col(t(mat))
  
  ga_entry = case_when(wk_entry == 1 ~ 16,
                       wk_entry == 2 ~ 17,
                       wk_entry == 3 ~ 18,
                       wk_entry == 4 ~ 19,
                       wk_entry == 5 ~ 20)
  
  #### Make the base dataset

  data <- tibble(n_sim, id, ga_entry, ptb_nocens, outcome, t_weeks) %>% 
    # Calculate time to event from ga_entry
    mutate(tte = t_weeks - ga_entry) %>% 
    # Exclude those individuals who experienced a delivery before trial entry
    filter(tte >= 0)
  
  return(data)
}


























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
