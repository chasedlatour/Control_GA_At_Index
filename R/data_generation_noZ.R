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
prob_enter <- c(0.2, 0.2, 0.2, 0.2, 0.2)
p_ptb0<- 0.2 # probability of delivery by 37 weeks gestation (i.e., preterm birth)
p_ptb1<- 0.1 # probability of delivery by 37 weeks gestation (i.e., preterm birth)
cenrate <- 0.1 # How much censoring
first_ga <- 15
weeks <- seq(15, 22, by=1)
p_loss <- 0.05
shape_loss <- 0.4
shape_del <- 20
shape_cens <- 0.2




generate_pregs <- function(sim_id, n, prob_enter, #p_loss, shape_loss,
                           p_ptb0, ptb_tte_wk_change, first_ga, weeks,
                           shape_del, shape_cens, scale_cens){
  
  set.seed(128467)
  
  # Record key identifiers
  n_sim <- sim_id   # identifier associated with the simulated dataset
  id <- 1:n         # person-level identifier
  
  
  # ### GENERATE MISCARRIAGE OUTCOMES STARTING AT THE FIRST WEEK OF GESTATION
  # p_event <- p_loss # probability of an event (loss)
  # losstau <- (21-first_ga)/52 # Maximum time to miscarriage, transform to year timescale
  # alpha <- shape_loss # weibull shape for tiem from first GA to pregnancy loss
  # lambda <- (-log(1-p_event))^(1/alpha)/losstau
  # h <- exp(log(lambda))
  # 
  # # Generate the time to miscarriage from first_ga
  # time_to_miscarriage <- rweibull(n, shape = alpha, scale = 1/h)
  
  ### GENERATE DELIVERY OUTCOMES STARTING AT THE FIRST WEEK OF GESTATION
  
  # Create parameters for preterm birth generation
  p_event <- p_ptb0 # probability of an event (preterm birth) 
  birthtau <- (37-first_ga)/52 #Maximum pregnancy duration in years
  alpha <- shape_del # weibull shape for time from 16 weeks of gestation to outcome 0.4
  lambda <- (-log(1-p_event))^(1/alpha)/birthtau   # weibull scale for time to PTB
  h <- exp(log(lambda)) 
  
  # Generate the time to delivery from first_ga
  time_to_event <- rweibull(n, shape=alpha, scale=1/h)        # time to delivery from first_ga
  
  # Create an indicator for PTB in a setting without censoring -- WE INCLUDE MISCARRIAGE
  # IN THE DEFINITION OF PTB FOR THIS EXAMPLE
  # ptb_nocens_t0 <- ifelse(pmin(time_to_event,time_to_miscarriage) > birthtau, 0, 1)
  # tte_ptb_nocens_t0 <- ifelse(ptb_nocens_t0 == 1, pmin(time_to_event, time_to_miscarriage), birthtau)
  ptb_nocens_t0 <- ifelse(time_to_event > birthtau, 0, 1)
  tte_ptb_nocens_t0 <- ifelse(ptb_nocens_t0 == 1, time_to_event, birthtau) # from first_ga
  
  
 
  ### GENERATE THEIR ACTUAL GESTATIONAL WEEK AT ENROLLMENT, independent of when they had an outcome
  
  mat <- rmultinom(n=n, size=1, prob_enter)
  ga_entry <- weeks[max.col(t(mat))]
  
  
  ### DETERMINE IF THE PERSON WOULD HAVE ENTERED THE TRIAL
  exclude <- ifelse((tte_ptb_nocens_t0*52) < ga_entry, 1, 0) # Transform tte_ptb_nocens_t0 to weeks timescale
  
  
  
  ## GENERATE THEIR POST-TREATMENT DELIVERY TIME
  # We assume that treatment delays the time to delivery by ptb_tte_wk_change weeks
  
  # First determine the time-to-delivery from ga_entry
  time_to_del <- ((first_ga/52) + time_to_event) - ga_entry/52 
      # # Calculate the gestational age at delivery and then Subtract the gestational age at trial entry
  # Determine the time time-to-delivery from ga_entry after treatment (i.e., ga_entry)
  time_to_event_t1 <- time_to_del + ptb_tte_wk_change/52
  
  # Create an indicator for PTB in a setting without censoring under treatment
  ptb_nocens_t1 <- ifelse(time_to_event_t1 > ((37-ga_entry)/52), 0, 1)
  tte_ptb_nocens_t1 <- ifelse(ptb_nocens_t1 == 1, time_to_event_t1, ((37-ga_entry)/52))
  
  
  ### GENERATE LTFU AT THE START OF ENROLLMENT FOR EACH PERSON, NOT 16 WEEKS GESTATION
  
  # Generate the time to censoring
  time_to_censor <- rweibull(n, shape=shape_cens, scale=scale_cens)
  
  ### NOW DETERMINE THEIR OBSERVED OUTCOME, APPLYING TIME TO CENSORING AFTER ENROLLMENT
  time_to_censor_ga <- (ga_entry - first_ga)/52 + time_to_censor
  
  
  
  # Observed event times and types based on ga at enrollment
  tte <- pmin(tte_ptb_nocens_t1, time_to_censor_ga, ((37-ga_entry)/52))
  t_weeks <- first_ga + (tte * 52)   # Gestational age at end of follow-up
  outcome <- ifelse(ptb_nocens_t1 == 1 & tte == tte_ptb_nocens_t1,
                    1,
                    ifelse(
                      ptb_nocens_t1 == 0 & tte == ((37-ga_entry)/52),
                      2,
                      0
                    ))
  
 
  #### Make the base dataset

  data <- tibble(n_sim, id, ga_entry, ptb_nocens_t0, ptb_nocens_t1, outcome, t_weeks) %>% 
    # Calculate time to event from ga_entry
    mutate(tte = t_weeks - ga_entry) %>% 
    # Exclude those individuals who experienced a delivery before trial entry
    filter(exclude == 0)
  
  return(data)
}


# test <- generate_pregs(sim_id = 1, n=20000, #prob_enter = rep(0.125, 8),
#                        prob_enter=c(0.25, 0.2, 0.2, 0.05, 0.05, 0.05, 0.1, 0.1),
#                        p_loss = 0.2, shape_loss = 0.33,
#                p_ptb = 0.13, first_ga = 15, weeks = seq(15,22, by=1),
#                shape_del = 20, # 0.4
#                shape_cens = 0.4, scale_cens = 10) # 12, 0.4
#                # shape_cens = 3, scale_cens = 0.7)

test <- generate_pregs(sim_id = 1, n=20000, prob_enter=c(0.05, 0.05, 0.05, 0.1, 0.1, 0.2, 0.2, 0.25),
                       p_ptb0 = 0.5, ptb_tte_wk_change = 1, first_ga = 15, weeks = seq(15,22, by=1),
                       shape_del = 10, shape_cens = 0.4, scale_cens = 10) # 12, 0.4





prop.table(table(test$outcome))
table(test$ga_entry)
mean(test$ptb_nocens_t1)
mean(test$ptb_nocens_t0)
prop.table(table(test$ga_entry, test$outcome))
prop.table(table(test$ga_entry[test$outcome == 0]))
# Look at outcome risk by GA at entry
test %>% 
  group_by(ga_entry) %>% 
  summarize(
    risk_t1 = mean(ptb_nocens_t1),
    risk_t0 = mean(ptb_nocens_t0)
  )






# Estimate the probability of preterm birth using an AJ estimator, not adjusting for gestational age at study entry

mod_noga <- summary(survfit(Surv(tte, factor(outcome)) ~ 1, data=test))
# Now get the cumulative incidence overall
prob <- tail(data.frame(r = mod_noga$pstate[,"1"]), 1)

# Now estimate risks within each GA stratum and standardize

mod <- summary(survfit(Surv(tte, factor(outcome)) ~ ga_entry, data=test))
strat_risk <- tibble(ga = mod$strata, r = mod$pstate[,"1"]) %>%
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
  summarize(risk = sum(std_risk),
            nonstandardized = mean(r))


# output side by side

mean(test$ptb_nocens_t1) # Probability without censoring (truth)
prob
standardized


















library(dplyr)
library(survival)
library(tibble)

# 1) pooled (unstratified) AJ — final cumulative incidence for cause = 1
mod_overall <- summary(survfit(Surv(tte, factor(outcome)) ~ 1, data = test))
pooled_final <- tail(mod_overall$pstate[ , "1"], 1)

# 2) stratified AJ: need to align rows correctly
mod_strat <- summary(survfit(Surv(tte, factor(outcome)) ~ ga_entry, data = test))

# mod_strat$strata is a named int vector: names are like "ga_entry=15", values are row counts
rows_per_stratum <- as.integer(mod_strat$strata)
stratum_names <- names(mod_strat$strata)

# Expand stratum labels to match the rows of mod_strat$time / pstate
stratum_for_row <- rep(stratum_names, times = rows_per_stratum)

df_strat_rows <- tibble(
  time = mod_strat$time,
  stratum = stratum_for_row,
  cif = mod_strat$pstate[ , "1"]   # cause-specific pstate column for cause=1
)

# Get final CIF per stratum
final_by_stratum <- df_strat_rows %>%
  group_by(stratum) %>%
  summarize(final_cif = last(cif), .groups = "drop") %>%
  mutate(ga_entry = as.integer(sub("ga_entry=", "", stratum))) %>%
  arrange(ga_entry)

# 3) sample weights for ga_entry (ensure type match)
props <- test %>%
  group_by(ga_entry) %>%
  summarize(prop = n()/nrow(test), .groups = "drop")

# join and compute direct-standardized final CIF
final_df <- final_by_stratum %>%
  left_join(props, by = "ga_entry")

# If any strata in final_df have NA prop (no observations), handle appropriately:
final_df %>% print(n = Inf)

standardized_final <- sum(final_df$final_cif * final_df$prop, na.rm = TRUE)

# 4) Report and compare
tibble(
  pooled_final = pooled_final,
  standardized_final = standardized_final,
  diff = pooled_final - standardized_final
)





















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





library(survival)
# pooled (unstratified) AJ
mod_overall <- summary(survfit(Surv(tte, factor(outcome)) ~ 1, data = test))
pooled_final <- tail(mod_overall$pstate[, "1"], 1)  # final pooled CIF for cause=1

# your standardized pipeline should have produced std_at_t0 (or similar)
# (if not, compute it with the function below)
pooled_final; standardized
pooled_final - standardized[1,1]




# stratified AJ
mod_strat <- summary(survfit(Surv(tte, factor(outcome)) ~ ga_entry, data = test))

# Build stratum-specific final CIFs (cause 1)
rows_per_stratum <- mod_strat$strata
strata_names <- names(mod_strat$strata)
labels <- rep(strata_names, times = rows_per_stratum)
df_mod <- data.frame(time = mod_strat$time, strata = mod_strat$strata, cif = mod_strat$pstate[, "1"])
library(dplyr)
final_by_stratum <- df_mod %>%
  group_by(strata) %>%
  summarize(final_cif = last(cif), n_rows = n()) %>%
  arrange(strata)
final_by_stratum

props <- test %>% count(ga_entry) %>% mutate(prop = n / sum(n))
props
# If you used equal weighting (mean(r)) vs sample weighting (sum(r*prop)), compute both:
strat_risk <- final_by_stratum %>%
  mutate(ga_entry = as.integer(sub("ga_entry=", "", strata))) %>%
  left_join(props, by = c("ga_entry" = "ga_entry")) %>%
  mutate(prop = ifelse(is.na(prop), 0, prop))

weighted_avg <- sum(strat_risk$final_cif * strat_risk$prop)
equal_avg    <- mean(strat_risk$final_cif)
weighted_avg; equal_avg

library(ggplot2)
library(tidyr)
library(zoo)

# Expand stratum CIFs to full grid and forward-fill (same approach as earlier)
df_mod <- df_mod %>% mutate(ga_entry = as.integer(sub("ga_entry=", "", strata)))
full_times <- sort(unique(df_mod$time))
strata_vals <- sort(unique(df_mod$ga_entry))
df_full <- tidyr::expand_grid(time = full_times, ga_entry = strata_vals) %>%
  left_join(df_mod %>% select(time, ga_entry, cif), by = c("time","ga_entry")) %>%
  arrange(ga_entry, time) %>%
  group_by(ga_entry) %>%
  mutate(cif = zoo::na.locf(c(ifelse(is.na(cif[1]), 0, 0), cif), na.rm = FALSE)[-1]) %>%
  ungroup()

# standardized curve (weights = sample props)
df_full <- df_full %>% left_join(props, by = c("ga_entry" = "ga_entry"))
std_curve <- df_full %>% mutate(w = prop*cif) %>% group_by(time) %>% summarize(std_cif = sum(w, na.rm = TRUE))

# pooled curve
pooled_df <- data.frame(time = mod_overall$time, pooled_cif = mod_overall$pstate[, "1"])

ggplot() +
  geom_step(data = df_full, aes(x = time, y = cif, group = factor(ga_entry)), alpha = 0.3) +
  geom_step(data = pooled_df, aes(x = time, y = pooled_cif), color = "blue", size = 1) +
  geom_step(data = std_curve, aes(x = time, y = std_cif), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Stratum CIFs (light), pooled AJ (blue), standardized (red dashed)",
       x = "time", y = "CIF (cause=1)")








ggplot() +
  # Stratum CIFs with color mapped to ga_entry
  geom_step(
    data = df_full,
    aes(x = time, y = cif, color = factor(ga_entry), group = factor(ga_entry)),
    alpha = 0.6
  ) +
  
  # Pooled AJ curve (fixed blue)
  geom_step(
    data = pooled_df,
    aes(x = time, y = pooled_cif),
    color = "blue",
    size = 1
  ) +
  
  # Standardized curve (fixed red dashed)
  geom_step(
    data = std_curve,
    aes(x = time, y = std_cif),
    color = "red",
    linetype = "dashed",
    size = 1
  ) +
  
  labs(
    title = "Cumulative Incidence of Outcome",
    subtitle = "Stratum CIFs by GA entry (colors), pooled AJ (blue), standardized (red dashed)",
    x = "Time",
    y = "CIF (cause = 1)",
    color = "GA Entry"
  ) +
  
  theme_minimal()
