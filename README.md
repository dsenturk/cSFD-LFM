# cSFD-LFM

## Contents of the cML-LFM folder

- cML-LFM_Tutorial.Rmd: A step-by-step implementation of cML-LFM and associated procedures including data generating process of the first simulation scenario and evaluation of model parameter estimation. Details of simulation design, estimation of model parameters and performance evaluation are described in 'Contrastive Latent Models for Structured Functional Data'.
- Simulation_dgp.R: Functions for generating contrastive multilevel functional data sets with latent components and score variances of Scenario 1 described in the simulation section.
- cML-LFM_Functions.R: Functions for conducting the proposed estimation algorithm described in Algorithm 1.

## Contents of the cMV-LFM folder

- cMV-LFM_Tutorial.Rmd: A step-by-step implementation of cMV-LFM and associated procedures including data generating process of the first simulation scenario and evaluation of model parameter estimation. Details of simulation design, estimation of model parameters and performance evaluation are described in 'Contrastive Latent Models for Structured Functional Data'.
- Simulation_dgp.R: Functions for generating contrastive multivariate functional data sets with latent components and score variances of Scenario 1 described in the simulation section.
- cMV-LFM_Functions.R: Functions for conducting the proposed estimation algorithm described in Algorithm 2.


## Introduction

The contents of this folder allow for the implementation of cSFD-LFM for contrastive functional data settings with multilevel or multivariate structure as proposed in "Contrastive Latent Models for Structured Functional Data". Users can simulate contrastive data pairs as described in the Simulation section and apply the proposed estimation algorithm to fit the proposed model frameworks. Detailed instructions on how to perform the aforementioned procedures are included in the two Tutorial.Rmd files.

## Requirements

- The included R programs require R 4.4.2 and the packages and files listed in cML-LFM_Tutorial.Rmd and cMV-LFM_Tutorial.Rmd.

## Installation

Load the R program files into the global environment and install the required packages using commands in cML-LFM_Tutorial.Rmd and cMV-LFM_Tutorial.Rmd.
