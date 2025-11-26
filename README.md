
# gaussian-process-regression
This repository contains the code related to a gaussian process regression approach to filling in Line-P data.

# Line P Seasonal Trend Fitting

In this example, we fit seasonal cycles and linear trends to Line P ocean timeseries data using the matlab function fitrgp.

The general idea is to obtain a mean seasonal cycle and trend from a timeseries as well as an associated uncertainty.


**Purpose:** Educational example for the Line P Hackathon (Dec 2025)

---

## Overview

This repository contains implementations (Python & MATLAB) for analyzing ocean timeseries data. 

The **Matlab** directory contains 
- LineP_SeasonalTrend_Fitting_Example.m. This Matlab file contains a working example where a mean trend and seasonal cycle are estimated to a timeseries. 
- LineP_DIC_UMOL_KG_timeseries_1990_2019_10m.nc: 10-m timesereis of DIC from Line P-P26. Input for LineP_SeasonalTrend_Fitting_Example.m

The **Python** directory contains: 
- fit_gpr_mauna_loa_gpss2025.ipynb. This file contains basic concepts of how a GPR would be adjusted to a timeseries that contains a positive trend and a 1 yr periodic cycle. Lines of code to download the data are included. This could be adapted to an ocean timeseries. 

Both Matlab and Python examples work a bit similarly: 
- **Fit a model:** constant + linear trend + annual seasonal cycle
- **Estimate uncertainties:** using the estimated covariance matrix 
- **Visualize results:** draw 1000 samples from a multivariate normal distribution of the fitted model coefficients

## Some resources
- The **intro lectures** from **[this](https://gpss.cc/gpss24/program)** Gaussian Process Summer School 2024 are really informative. 
- The **Labs** from the 2025 summer school can be run online on a browser **[Lab 1: Gaussian Process Regression](https://colab.research.google.com/github/gpschool/gpss25/blob/gh-pages/labs/lab_1.ipynb#scrollTo=iXI-p5IKwzJv)**
- Rasmussen and Williams in their book **[Gaussian Processes for Machine Learning](https://gaussianprocess.org/gpml/chapters/RW.pdf)** tell you all you didn't know that you didn't know about GPR
