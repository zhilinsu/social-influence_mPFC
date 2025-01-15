# ============================================================================.
# Info #### 
# ============================================================================.
# Preference-uncertainty (KU) model for multiple subjects. 
# Hierarchical structure with uninformative priors.
# 
# Ref: Moutoussis et al., 2016, *PLoS Comput. Biol.*  
# 
# Zhilin Su
# zhilinsu1312@gmail.com / z.su.1@pgr.bham.ac.uk 
#
# Project: Zhilin's PhD project 2 - mPFC and social influence.   
# Last revision: 15 Jan 2025, Zhilin Su. 

# ============================================================================.
# Preparation #### 
# ============================================================================.
# Clear all variables 
rm(list = ls())

# Load libraries
library(rstan)
library(ggplot2)
library(dplyr)
library(abind)

# Set seed 
set.seed(1213)

# ============================================================================.
# Set Populations & Initialisation #### 
# ============================================================================.
# 0 for combining all groups together, 
# 1 for healthy controls, 3 for mPFC lesions, and 4 for lesion controls. 
pop_vec <- c(1, 3, 4)

for (pop_value in pop_vec) {
	
	pop <- pop_value
	
	load("data/sd_data_mpfc.RData")
	
	if (pop == 0) {
		pt_1 <- sd_hc
		pt_2 <- sd_mpfc
		pt_3 <- sd_lc
		other_order_1 <- k_hc 
		other_order_2 <- k_mpfc 
		other_order_3 <- k_lc
		incl_id_1 <- readRDS("data/incl_id_hc.rds")
		incl_id_2 <- readRDS("data/incl_id_mpfc.rds")
		incl_id_3 <- readRDS("data/incl_id_lc.rds")
		
		pt <- abind(pt_1, pt_2, pt_3, along = 3)
		other_order <- rbind(other_order_1, other_order_2, other_order_3)
		incl_id <- c(incl_id_1, incl_id_2, incl_id_3)
		
	} else if (pop == 1) {
		pt <- sd_hc 
		other_order <- k_hc
		incl_id <- readRDS("data/incl_id_hc.rds")
		
	} else if (pop == 3) {
		pt <- sd_mpfc
		other_order <- k_mpfc 
		incl_id <- readRDS("data/incl_id_mpfc.rds")
		
	} else if (pop == 4) {
		pt <- sd_lc
		other_order <- k_lc 
		incl_id <- readRDS("data/incl_id_lc.rds")
	}
	
	rm(sd_hc, sd_mpfc, sd_lc, k_hc, k_mpfc, k_lc)
	
	# Convert data type 
	other_order <- other_order %>% 
		as_tibble() %>%
		mutate(id = as.numeric(id), 
					 other1_k = as.numeric(other1_k), 
					 other2_k = as.numeric(other2_k), 
					 other1_pref = as.factor(other1_pref), 
					 other2_pref = as.factor(other2_pref))
	
	# Include participants 
	incl_condition <- pt[1, "id", ] %in% incl_id 
	pt <- pt[,, incl_condition]
	other_order <- other_order %>% 
		filter(id %in% incl_id)
	
	# Extract metadata
	n_trials <- dim(pt)[1] # number of trials 
	n_pt <- dim(pt)[3] # number of participants
	
# ============================================================================.
# Reshape data frames #### 
# ============================================================================.
	id_vec <- vector("numeric", length = n_pt)
	id_vec <- pt[1, "id", 1:n_pt]
	names(id_vec) <- NULL
	
	# Initialise a 3D array 
	pt_new <- array(0, dim = dim(pt))
	col_names <- c("id", "s_dec_1", "s_rs_1", "s_rl_1", "s_dl_1", 
								 "so_dec_i", "o_rs_i", "o_rl_i", "o_dl_i", "o_dec_i", 
								 "s_dec_i", "s_rs_i", "s_rl_i", "s_dl_i", 
								 "so_dec_p", "o_rs_p", "o_rl_p", "o_dl_p", "o_dec_p", 
								 "s_dec_p", "s_rs_p", "s_rl_p", "s_dl_p")
	colnames(pt_new) <- col_names
	
	# Reshape the data 
	for (i in 1:n_pt) {
		other1_pref <- other_order$other1_pref[i]
		
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
# Set Stan Parameters #### 
# ============================================================================.
	rstan_options(auto_write = TRUE)
	options(mc.cores = 4)
	
	model_file <- "scripts/ku0-model_self.stan"
	
	n_iter     <- 4000
	n_chains   <- 4 
	n_warmup   <- floor(n_iter / 2)
	n_thin     <- 1
	n_seed     <- 1213 # set the seed
	
# ============================================================================.
# Main Loop for Running Stan #### 
# ============================================================================.
	for (s_block in c("1", "i", "p")) { # self blocks 
		# 1 - Construct data
		## Initialise the arrays 
		self_decision <- array(0, dim = c(n_pt, n_trials)) 
		reward_sooner <- array(0, dim = c(n_pt, n_trials))
		reward_later <- array(0, dim = c(n_pt, n_trials))
		delay_later <- array(0, dim = c(n_pt, n_trials))
		
		for (i in 1:n_pt) {
			self_decision[i, ] <- pt_new[, paste("s_dec_", s_block, sep = ""), i] 
			reward_sooner[i, ] <- pt_new[, paste("s_rs_", s_block, sep = ""), i]
			reward_later[i, ] <- pt_new[, paste("s_rl_", s_block, sep = ""), i]
			delay_later[i, ] <- pt_new[, paste("s_dl_", s_block, sep = ""), i] 
		}
		
		data_list <-  list(n_trials   = n_trials,
											 n_subjects = n_pt,
											 s_dec      = self_decision,
											 s_rs       = reward_sooner,
											 s_rl       = reward_later,
											 s_dl       = delay_later)
			
		# 2 - Run Stan 
		cat("Estimating", model_file, "model... \n")
		cat("Estimating for the Block", s_block, "... \n")
		start_time <- Sys.time()
		print(start_time)
		cat("Calling", n_chains, "simulations in Stan... \n")
		
		fit_kt <- stan(model_file, 
									 data   = data_list,
									 chains = n_chains, 
									 iter   = n_iter, 
									 warmup = n_warmup, 
									 thin   = n_thin,
									 init   = "random",
									 seed   = n_seed)
		
		cat("Finishing", model_file, "model simulation ... \n")
		cat("Finishing the estimation for the Block", s_block, "... \n")
		end_time <- Sys.time()
		print(end_time)  
		cat("It took", as.character.Date(end_time - start_time), "\n")
		
		# 3 - Save stanfit objects 
		if (pop == 0) {
			file_name <- paste("stanfit_ku0_all_self_", s_block, ".rds", sep = "")
		} else if (pop == 1) {
			file_name <- paste("stanfit_ku0_hc_self_", s_block, ".rds", sep = "")
		} else if (pop == 3) {
			file_name <- paste("stanfit_ku0_mpfc_self_", s_block, ".rds", sep = "")
		} else if (pop == 4) {
			file_name <- paste("stanfit_ku0_lc_self_", s_block, ".rds", sep = "")
		}
		
		saveRDS(fit_kt, paste("data/stanfit/", file_name, sep = "")) 
		
		# 4 - Clear variables 
		rm(data_list, fit_kt)
	}
}