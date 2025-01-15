Code and data for the delegated inter-temporal choice task, accompanying the manuscript:

Zhilin Su\*, Mona M. Garvert, Lei Zhang, Todd A. Vogel, Jo Cutler, Masud Husain, Sanjay Manohar, & Patricia L. Lockwood\*. (*under revision*). **Dorsomedial and ventromedial prefrontal cortex lesions differentially impact social influence and temporal impulsivity**.

This repository contains:

```
root
  ├─ data     # behavioural data and plots 
       ├─ plots 
       ├─ stanfit # stanfit objects
  ├─ lesion   # lesion analysis
       ├─ all_lesion_patients 
       ├─ mPFC_patients_only
  ├─ scripts  # R and Stan code to run the analyses and produce figures
       ├─ helper_functions 
       ├─ stan_model 
```

# Installation

All behavioural analyses were performed using R (v4.2.1) in RStudio (v2023.06.2). Installation guides can be found at <https://posit.co/download/rstudio-desktop/>.

Model fitting was performed using Stan (v2.32) and the RStan (v2.21.7) package in RStudio. Installation guides can be found at <https://mc-stan.org/users/interfaces/>.

All these installations should only take a few minutes to complete.

# Preparation 

**0 - Data**

Data for the delegated inter-temporal choice task collected through MATLAB have been extracted and aggregated into `sd_data_mpfc.RData`, for the subsequent analysis and modelling using R and Stan. 

Other self-reported data were stored in `demographics.csv` and `questionnaires.csv`.

# Model fitting 

**1 - Bayesian modelling with Stan**

Run the R scripts `run_*1_model_*2.R` to call the corresponding Stan scripts to perform Bayesian modelling.

**2 - Parameter recovery**

**3 - Posterior predictive checks**

# Behavioural analysis 

# Lesion analysis 