# README
Valdrich Fernandes, Perry de Louw, Coen Ritsema, Ruud Bartholomeus

## Introduction

This repository consists of the code used to analyse and predict the
transient response to Managed Aquifer Recharge (MAR). The results of
this code is under review, under a publication named **“Predicting
Groundwater Storage from Seasonal MAR: Insights from Machine Learning
and Interpretable AI Techniques”**.

## Abstract

Managed Aquifer Recharge (MAR) is widely applied to enhance groundwater
storage and promote the sustainable use of this essential resource.
Techniques such as Machine Learning (ML) as surrogate for
computationally intensive numerical models, are increasingly considered
for identifying suitable locations for MAR. While ML has been
demonstrated to be suitable for steady-state simulations, its
application for transient modelling is much more challenging. However,
understanding how water recharged during wet periods is retained and
remains available during dry periods is critical, highlighting the
importance of transient responses. In this study, we therefore employed
ML to mimic the transient effect of MAR by decomposing the time series
of groundwater storage when MAR stopped into two components—the
MAR-response and the decay coefficient—following an exponential decay
curve. This decomposition provides a simplified yet effective
representation of groundwater storage changes over time. Using U-Net and
XGBoost, we demonstrated that ML can accurately capture these dynamics
for the Baakse Beek catchment, located in the drought sensitive sandy
soils of the Netherlands. The ML models achieved an R2 of more than 0.85
in predicting the two components of the long-term dynamics of the stored
water. Additionally, using explainable AI techniques, specifically
SHapley Aditive exPlanation values, we identified the site management
decisions and the surface water network properties near the recharge
site most significantly impact the effectiveness of the MAR in the
region. This focus on model interpretability ensures transparency in the
predictions, improving on the models’ generalizability and fosters trust
among hydrologists and stakeholders.
