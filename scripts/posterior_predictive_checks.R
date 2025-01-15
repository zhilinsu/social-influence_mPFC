# ============================================================================.
# Info ####
# ============================================================================.
# The script is to perform posterior predictive checks of the KU model. 
# 
# Zhilin Su
# zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
#
# Project: Zhilin's second PhD project - mPFC and social influence
# Last revision: 15 January 2025, Zhilin Su 

rm(list = ls())

library(rstan)
library(tidyverse)
library(ggtext)
library(lme4)
library(parameters)
library(kableExtra)
library(ggsignif)
library(BayesFactor)
source("scripts/helper_functions/helper_functions.R")

set.seed(1213)

# 1 for HC impulsive; 2 for mPFC impulsive; 3 for LC impulsive; 
# 4 for HC patient; 5 for mPFC patient; 6 for LC patient; 
acc_result_list <- vector(mode = "list", length = 6)

# 1 for healthy controls, 3 for mPFC patients, 4 for lesion controls 
pop <- 4

# ============================================================================.
# Preparation #### 
# ============================================================================.
load("data/sd_data_mpfc.RData")

if (pop == 1) {
	pt <- sd_hc
	other_k <- k_hc 
	incl_id <- readRDS("data/incl_id_hc.rds")
} else if (pop == 3) {
	pt <- sd_mpfc
	other_k <- k_mpfc 
	incl_id <- readRDS("data/incl_id_mpfc.rds")
} else if (pop == 4) { 
	pt <- sd_lc
	other_k <- k_lc 
	incl_id <- readRDS("data/incl_id_lc.rds")
	}

# Pre-process other_order 
other_k <- other_k %>% 
	as_tibble() %>%
	mutate(id = as.numeric(id), 
				 other1_k = as.numeric(other1_k), 
				 other2_k = as.numeric(other2_k), 
				 other1_pref = as.factor(other1_pref), 
				 other2_pref = as.factor(other2_pref)) %>%
	mutate(other_i_k = if_else(other1_pref == "i", other1_k, other2_k),
				 other_p_k = if_else(other1_pref == "p", other1_k, other2_k))

# Include participants 
incl_condition <- pt[1, "id", ] %in% incl_id 
pt <- pt[,, incl_condition]
other_k <- other_k %>% 
	filter(id %in% incl_id)

# Extract metadata
sz <- dim(pt) # sz = [# of trials, # of experiment parameters, # of participants]
n_trials <- sz[1] # number of trials 
n_pt <- sz[3] # number of participants

# Reshape data 
id_vec <- vector("numeric", length = n_pt)
id_vec <- pt[1, "id", 1:n_pt]
names(id_vec) <- NULL

if (pop == 1) {
	id_vec_hc <- id_vec
} else if (pop == 3) {
	id_vec_mpfc <- id_vec
} else if (pop == 4) { 
	id_vec_lc <- id_vec
}

pt_new <- array(0, dim = sz)
col_names <- c("id", "s_dec_1", "s_rs_1", "s_rl_1", "s_dl_1", 
							 "so_dec_i", "o_rs_i", "o_rl_i", "o_dl_i", "o_dec_i", 
							 "s_dec_i", "s_rs_i", "s_rl_i", "s_dl_i", 
							 "so_dec_p", "o_rs_p", "o_rl_p", "o_dl_p", "o_dec_p", 
							 "s_dec_p", "s_rs_p", "s_rl_p", "s_dl_p")
colnames(pt_new) <- col_names

# Reshape data 
for (i in 1:n_pt) {
	
	other1_pref <- other_k$other1_pref[i]
	
	if (other1_pref == "i") {
		pt_new[, 1:5, i] <- pt[, 1:5, i] # from id to s_dl_1
		pt_new[, 6:14, i] <- pt[, 6:14, i] # from so_dec_i to s_dl_i 
		pt_new[, 15:23, i] <- pt[, 15:23, i] # from so_dec_p to s_dl_p 
	} else if (other1_pref == "p") {
		pt_new[, 1:5, i] <- pt[, 1:5, i] # from id to s_dl_1 
		pt_new[, 6:14, i] <- pt[, 15:23, i] # from so_dec_i to s_dl_i 
		pt_new[, 15:23, i] <- pt[, 6:14, i] # from so_dec_p to s_dl_p 
	}
} 

# ============================================================================.
# Set Stan parameters ####
# ============================================================================.
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_iter     <- 4000
n_chains   <- 4 
n_warmup   <- floor(n_iter / 2)
n_thin     <- 1
n_seed     <- 1213

# ============================================================================.
# Data simulation #### 
# ============================================================================.
## Prepare experimental parameters ####
# Impulsive block 
i_reward_sooner_sim <- array(0, dim = c(n_pt, n_trials))
i_reward_later_sim <- array(0, dim = c(n_pt, n_trials))
i_delay_later_sim <- array(0, dim = c(n_pt, n_trials))

# Patient block 
p_reward_sooner_sim <- array(0, dim = c(n_pt, n_trials))
p_reward_later_sim <- array(0, dim = c(n_pt, n_trials))
p_delay_later_sim <- array(0, dim = c(n_pt, n_trials))

for (i in 1:n_pt) {
	i_reward_sooner_sim[i, ] <- pt_new[, "o_rs_i", i]
	i_reward_later_sim[i, ] <- pt_new[, "o_rl_i", i]
	i_delay_later_sim[i, ] <- pt_new[, "o_dl_i", i]
	
	p_reward_sooner_sim[i, ] <- pt_new[, "o_rs_p", i]
	p_reward_later_sim[i, ] <- pt_new[, "o_rl_p", i]
	p_delay_later_sim[i, ] <- pt_new[, "o_dl_p", i]
}

## Get empirical posterior values 
km_i <- array(0, n_pt)
ku_i <- array(0, n_pt)

km_p <- array(0, n_pt)
ku_p <- array(0, n_pt)

if (pop == 1) {
	ku_model_other_i_ya <- 
		readRDS("data/stanfit/stanfit_ku0_hc_other_i.rds")
	ku_model_other_p_ya <- 
		readRDS("data/stanfit/stanfit_ku0_hc_other_p.rds")
	
	km_i <- colMeans(rstan::extract(ku_model_other_i_ya)$km)
	ku_i <- colMeans(rstan::extract(ku_model_other_i_ya)$ku)
	
	km_p <- colMeans(rstan::extract(ku_model_other_p_ya)$km)
	ku_p <- colMeans(rstan::extract(ku_model_other_p_ya)$ku)
	
} else if (pop == 3) {
	ku_model_other_i_oa <- 
		readRDS("data/stanfit/stanfit_ku0_mpfc_other_i.rds")
	ku_model_other_p_oa <- 
		readRDS("data/stanfit/stanfit_ku0_mpfc_other_p.rds")
	
	km_i <- colMeans(rstan::extract(ku_model_other_i_oa)$km)
	ku_i <- colMeans(rstan::extract(ku_model_other_i_oa)$ku)
	
	km_p <- colMeans(rstan::extract(ku_model_other_p_oa)$km)
	ku_p <- colMeans(rstan::extract(ku_model_other_p_oa)$ku)
	
} else if (pop == 4) {
	ku_model_other_i_oa <- 
		readRDS("data/stanfit/stanfit_ku0_lc_other_i.rds")
	ku_model_other_p_oa <- 
		readRDS("data/stanfit/stanfit_ku0_lc_other_p.rds")
	
	km_i <- colMeans(rstan::extract(ku_model_other_i_oa)$km)
	ku_i <- colMeans(rstan::extract(ku_model_other_i_oa)$ku)
	
	km_p <- colMeans(rstan::extract(ku_model_other_p_oa)$km)
	ku_p <- colMeans(rstan::extract(ku_model_other_p_oa)$ku)
}

simulated_model_file <- 
	"scripts/ku0-model_other_simulated.stan"

i_sim_data_list <-  list(n_trials   = n_trials,
												 n_subjects = n_pt,
												 o_rs       = i_reward_sooner_sim,
												 o_rl       = i_reward_later_sim,
												 o_dl       = i_delay_later_sim,
												 km         = km_i, 
												 ku         = ku_i)

p_sim_data_list <-  list(n_trials   = n_trials,
												 n_subjects = n_pt,
												 o_rs       = p_reward_sooner_sim,
												 o_rl       = p_reward_later_sim,
												 o_dl       = p_delay_later_sim,
												 km         = km_p, 
												 ku         = ku_p)

if (pop == 1) {
	saveRDS(i_sim_data_list, "data/i_sim_data_list_hc.RDS")
	saveRDS(p_sim_data_list, "data/p_sim_data_list_hc.RDS")
} else if (pop == 3) {
	saveRDS(i_sim_data_list, "data/i_sim_data_list_mpfc.RDS")
	saveRDS(p_sim_data_list, "data/p_sim_data_list_mpfc.RDS")
} else if (pop == 4) {
	saveRDS(i_sim_data_list, "data/i_sim_data_list_lc.RDS")
	saveRDS(p_sim_data_list, "data/p_sim_data_list_lc.RDS")
}

i_sim_data <- stan(simulated_model_file,
									 data = i_sim_data_list,
									 chains = n_chains,
									 iter   = n_iter,
									 warmup = n_warmup,
									 thin   = n_thin,
									 init   = "random", 
									 seed   = n_seed, 
									 algorithm = "Fixed_param")

p_sim_data <- stan(simulated_model_file,
									 data = p_sim_data_list,
									 chains = n_chains,
									 iter   = n_iter,
									 warmup = n_warmup,
									 thin   = n_thin,
									 init   = "random",
									 seed   = n_seed,
									 algorithm = "Fixed_param")

if (pop == 1) {
	saveRDS(i_sim_data, "data/i_sim_data_hc.RDS")
	saveRDS(p_sim_data, "data/p_sim_data_hc.RDS")
} else if (pop == 3) {
	saveRDS(i_sim_data, "data/i_sim_data_mpfc.RDS")
	saveRDS(p_sim_data, "data/p_sim_data_mpfc.RDS")
} else if (pop == 4) {
	saveRDS(i_sim_data, "data/i_sim_data_lc.RDS")
	saveRDS(p_sim_data, "data/p_sim_data_lc.RDS")
}

# ============================================================================.
# Accuracy calculation #### 
# ============================================================================.
# rstan::extract(i_sim_data)$s_dec_sim[a, b, c] 
# a: iteration; b: participant; c: trial

value_later <- array(0, dim = c(n_pt, n_trials))
model_choice <- array(0, dim = c(n_pt, n_trials)) 
acc_result <- array(0, n_pt)

## For impulsive #### 
for (j in 1:n_pt) {
	k_value <- 10^(as.numeric(other_k[j, "other_i_k"]))
	value_later[j, ] <- i_reward_later_sim[j, ] / (1 + k_value * i_delay_later_sim[j, ])
	model_choice[j, ] <- (value_later[j, ] > i_reward_sooner_sim[j, ]) + 1
	
	acc_vec <- array(0, 8000)
	s_dec_sim_matrix <- rstan::extract(i_sim_data)$s_dec_sim[, j, ]
	
	for (i in 1:8000) {
		sim_dec <- s_dec_sim_matrix[i, ]
		acc_vec[i] <- mean(sim_dec == model_choice[j, ])
		
	} # iteration loop 
	acc_result[j] = mean(acc_vec)
} # participant loop 

if (pop == 1) {
	acc_result_list[[1]] <- acc_result 
} else if (pop == 3) {
	acc_result_list[[2]] <- acc_result 
} else if (pop == 4) {
	acc_result_list[[3]] <- acc_result
}

## For patient ####
for (j in 1:n_pt) {
	k_value <- 10^(as.numeric(other_k[j, "other_p_k"]))
	value_later[j, ] <- p_reward_later_sim[j, ] / (1 + k_value * p_delay_later_sim[j, ])
	model_choice[j, ] <- (value_later[j, ] > p_reward_sooner_sim[j, ]) + 1
	
	acc_vec <- array(0, 8000)
	s_dec_sim_matrix <- rstan::extract(p_sim_data)$s_dec_sim[, j, ]
	
	for (i in 1:8000) {
		sim_dec <- s_dec_sim_matrix[i, ]
		acc_vec[i] <- mean(sim_dec == model_choice[j, ])
		
	} # iteration loop 
	acc_result[j] = mean(acc_vec)
} # participant loop 

if (pop == 1) {
	acc_result_list[[4]] <- acc_result 
} else if (pop == 3) {
	acc_result_list[[5]] <- acc_result 
} else if (pop == 4) {
	acc_result_list[[6]] <- acc_result
}

saveRDS(acc_result_list, "data/acc_result_list.RDS")

# ============================================================================.
# Analyses #### 
# ============================================================================.
# Convert list to tibble 
sim_acc_result <- bind_rows(
	tibble(group = "1", preference = "i", accuracy = acc_result_list[[1]], 
				 id = id_vec_hc),
	tibble(group = "1", preference = "p", accuracy = acc_result_list[[4]],
				 id = id_vec_hc),
	tibble(group = "3", preference = "i", accuracy = acc_result_list[[2]],
				 id = id_vec_mpfc),
	tibble(group = "3", preference = "p", accuracy = acc_result_list[[5]],
				 id = id_vec_mpfc),
	tibble(group = "4", preference = "i", accuracy = acc_result_list[[3]],
				 id = id_vec_lc),
	tibble(group = "4", preference = "p", accuracy = acc_result_list[[6]],
				 id = id_vec_lc)
) %>% 
	mutate(group = factor(group, 
												levels = c("3", "1", "4"), 
												labels = c("mPFC", "HC", "LC")), 
				 preference = factor(preference), 
				 accuracy = as.numeric(accuracy)) %>%
	mutate(id = factor(id)) %>%
	mutate(accuracy = accuracy * 100) %>%
	mutate(condition = interaction(group, preference, sep = ":"))

saveRDS(sim_acc_result, "data/sim_acc_result.RDS")

sim_acc_result <- readRDS("data/sim_acc_result.RDS")

# LMM 
sim_acc_result <- sim_acc_result %>%
	mutate(group = factor(group, 
												levels = c("HC", "mPFC", "LC")))

contrasts(sim_acc_result$preference) <- c(-1, 1)

lmer_model <- lmer(accuracy ~ 1 + group * preference + (1|id), 
									 data = sim_acc_result)

acc_lmer_table_full <- model_parameters(lmer_model)

acc_lmer_table <- acc_lmer_table_full %>%
	filter(Effects == "fixed") %>% 
	dplyr::select(!c(df_error, Effects, Group, CI)) %>%
	mutate(Parameter = str_replace(Parameter, pattern = "groupHC", replacement = "Group (HC)")) %>%
	mutate(Parameter = str_replace(Parameter, pattern = "groupLC", replacement = "Group (LC)")) %>%
	mutate(Parameter = str_replace(Parameter, pattern = "other_pref1", replacement = "Others (patient)"))

acc_lmer_table

# ============================================================================.
# Visualisation #### 
# ============================================================================.
vars <- "accuracy"

plt_data <- sim_acc_result %>% 
	group_by(group, preference, condition) %>% 
	filter(!is.na(get(vars))) %>%
	summarise(n = n(), 
						mean = mean(get(vars), na.rm = TRUE), 
						sd = sd(get(vars), na.rm = TRUE), 
						se = sd / sqrt(n)) %>% 
	ungroup()

plt <- ggplot(plt_data, aes(x = group, y = mean, group = condition)) +
	geom_jitter(data = sim_acc_result,
							aes(x = group, y = get(vars), 
									colour = preference, 
									alpha = preference), 
							position = position_jitterdodge(dodge.width = 0.75, 
																							jitter.height = 0, 
																							jitter.width = 0.3), 
							size = 1, stroke = 0.5,
							show.legend = FALSE) + 
	geom_errorbar(aes(ymin = mean - se, 
										ymax = mean + se), 
								position = position_dodge(0.75), 
								colour = "black", width = 0, size = 0.25) + 
	geom_jitter(aes(y = mean, 
									fill = preference), 
							position = position_dodge(0.75), 
							colour = "black", shape = 21, stroke = 0.25, size = 2) +
	scale_x_discrete(labels = c("Healthy\ncontrols", 
															"mPFC\nlesions", 
															"Lesion\ncontrols")) +
	coord_cartesian(ylim = c(50, 105), expand = FALSE) +  
	labs(x = NULL, y = "Simulated learning accuracy (%)") +
	scale_fill_manual(name = "other_pref", 
										limits = c("i", "p"), 
										values = c("#5491cb", "#fcd466"), 
										labels = c("Impulsive", "Patient")) +
	scale_colour_manual(name = "other_pref",
											limits = c("i", "p"),
											values = c("#abdeff", "#ffeacb")) +
	scale_alpha_manual(name = "other_pref",
										 limits = c("i", "p"),
										 values = c(0.3, 0.3)) + 
	theme_classic() +
	theme(text = element_text(family = "Arial"),
				legend.title = element_blank(),
				legend.key = element_rect(linewidth = c(2, 2), colour = NA),
				axis.text.x = element_text(size = 10),
				axis.text.y = element_markdown(size = 10), 
				axis.title.x = element_markdown(size = 10), 
				axis.title.y = element_markdown(size = 12, lineheight = 1.2), 
				legend.text = element_text(size = 8, margin = margin(l = 0)),
				legend.key.height = unit(0.3, "cm"),
				legend.key.width = unit(0.3, "cm"),
				legend.key.spacing.y = unit(0.1, "cm"),
				legend.justification = c(0, 0.5)) +
	annotate("point", x = 0.9, y = 80.8, colour = "red", alpha = 0.75, size = 1) +
	annotate("point", x = 1.27, y = 85.1, colour = "red", alpha = 0.75, size = 1) +
	annotate("point", x = 1.9, y = 79.4, colour = "red", alpha = 0.75, size = 1) +
	annotate("point", x = 2.27, y = 80.3, colour = "red", alpha = 0.75, size = 1) +
	annotate("point", x = 2.9, y = 78.0, colour = "red", alpha = 0.75, size = 1) +
	annotate("point", x = 3.27, y = 80.3, colour = "red", alpha = 0.75, size = 1) +
	geom_signif(comparisons = list(c("HC", "mPFC")),
							y_position = 96.5, tip_length = 0.01, vjust = 0.5, size = 0.25,
							textsize = 3, family = "Arial",
							annotation = "***") + 
	geom_signif(comparisons = list(c("HC", "LC")),
							y_position = 92.5, tip_length = 0.01, vjust = 0.5, size = 0.25,
							textsize = 3, family = "Arial",
							annotation = "***")

ggsave("data/plots/simulated_acc_ku.png", plt,
			 height = 4, width = 4 * 0.85, dpi = 1200)
 
