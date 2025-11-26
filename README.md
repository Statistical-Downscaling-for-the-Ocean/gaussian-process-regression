
# gaussian-process-regression
This repository contains the code related to a gaussian process regression approach to filling in Line-P data.

# Line P Seasonal Trend Fitting

In this example, we fit seasonal cycles and linear trends to Line P ocean timeseries data using the matlab function fitrgp.

The general idea is to obtain a mean seasonal cycle and trend from a timeseries as well as an associated uncertainty.


**Purpose:** Educational example for the Line P Hackathon (Dec 2025)

---

## Overview

This repository contains implementations (Python & MATLAB) for analyzing ocean timeseries data:

- **Fit a model:** constant + linear trend + annual seasonal cycle
- **Estimate uncertainties:** using the estimated covariance matrix 
- **Visualize results:** draw 1000 samples from a multivariate normal distribution of the fitted model coefficients

## Some resources
- The **intro lectures** from **[this](https://gpss.cc/gpss24/program)** Gaussian Process Summer School 2024 are really informative. 
- The **Labs** from the 2025 summer school can be run online on a browser **[Lab 1: Gaussian Process Regression](https://colab.research.google.com/github/gpschool/gpss25/blob/gh-pages/labs/lab_1.ipynb#scrollTo=iXI-p5IKwzJv)**
- Rasmussen and Williams in their book **[Gaussian Processes for Machine Learning](https://gaussianprocess.org/gpml/chapters/RW.pdf)** tell you all you didn't know that you didn't know about GPR
