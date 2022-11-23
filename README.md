## Code and data associated with: 
# Quantifying the impact of immune history and variant on SARS-CoV-2 viral kinetics and infection rebound: a retrospective cohort study
Published at: https://doi.org/10.7554/eLife.81849

__James A. Hay, Stephen M. Kissler, Joseph R. Fauver, Christina Mack, Caroline G. Tai, Radhika M. Samant, Sarah Connolly, Deverick J. Anderson, Gaurav Khullar, Matthew MacKay, Miral Patel, Shannan Kelly, April Manhertz, Isaac Eiter, Daisy Salgado, Tim Baker, Ben Howard, Joel T. Dudley, Christopher E. Mason, Manoj Nair, Yaoxing Huang, David D. Ho, Nathan D. Grubaugh, Yonatan H. Grad__

For issues with code, please contact: 

jhay@hsph.harvard.edu

skissler@hsph.harvard.edu

ygrad@hsph.harvard.edu


All code is available under the GNU General Public License, version 3, included in this repository under the following terms: 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Summary
Code and data to reproduce all of the analyses in the paper linked here: [https://doi.org/10.1101/2022.01.13.22269257](https://doi.org/10.1101/2022.01.13.22269257).

The code is broadly split into two parts: 1) code to clean the data and generate the logistic regression fits using the `brms` package, and 2) code to run the piece-wise linear regression model of viral kinetics using `RStan`.

<<NOTE>> Due to the sensitive nature of these data, the raw files provided here have been stripped of private, potentially identifiable information, including all dates, age and role. Thus, some of the analyses and figures are not directly reproducible from the raw data files and the provided R scripts will not run without error. However, all of the data to reproduce viral kinetics as a function of immune history, variant and time-since-detection are included.

# Setup
Key packages to install are [`RStan`](https://mc-stan.org/users/interfaces/rstan) and [`brms`](https://cran.r-project.org/web/packages/brms/index.html). Both require an installation of Stan and a working C++ toolchain, installation of which depends on your operating system (see [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)).

Many of the analyses are run on the [Harvard FASRC HPC](https://www.rc.fas.harvard.edu/) as batch jobs. These are any of the scripts in the `cluster` folder. It is possible to tweak these to run locally, but it may require a bit of fiddling. In short, we run these analyses by closing the Git repo to a home working directory on the cluster and sourcing scripts/using batch files to submit jobs through the command line.

Otherwise, standard R packages are required including:
``` r
c("tidyverse","lubridate","ggbeeswarm","data.table","patchwork","mice","zoo",
"splines","brms","mgcv","ROCR","tidybayes","RStan")
```

# Data cleaning
Provided in the `data` folder are two versions of the raw data:
1. `ct_data_cleaned.csv`: the raw data prior to pre-processing for regression analyses. Note that there are a number of variables which are used only in these scripts.
2. `data_for_regressions.RData`: the cleaned data used for all of the regression analyses and viral kinetics model. This has some additional variables, transformations, filters and categorizations to facilitate the model fitting.

# Raw data analyses
`ct_data_cleaned.csv` is used to generate raw counts and summary statistics using the following scripts:
1. `prep_data_for_regression.R`: adds variables, transformations and filters to the raw data to facilitate model fitting. Also produces Fig 1A and Fig 2C and generates the symptomatic data statistics.
2. `raw_data.R`: generates the summary statistics for the results section and Fig S1.
3. `rebounds.R`: generates Fig 1B. Table 1 and the summary statistics on rebound infections. Also generates Figure 1.

# Logistic regression analyses
All of these analyses are run on the Harvard FASRC HPC. For each analysis, there are two scripts:
1. `cluster_submit_XX.sh`: this is the batch file used to submit a list of jobs.
2. `brm_regression_cluster_XX.R`: the R script sourced by the cluster submit file. Sourcing this file once will run a single iteration of the model fitting. You can modify this file to run locally by changing the working directory to this Git repository, and by setting `i` to the desired model run.

Outputs from these scripts are saved in the `outputs` folder. We have only included the RStan objects for a subset of the base regression models, immune history regression models and boost/titer group models due to size constraints. All of the remaining RStan objects, including the k-folds cross-validation dataset will be provided if requested.

# Viral kinetics model
All code for this part of the analyses is in the `code` folder.
