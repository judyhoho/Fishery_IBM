# Code for: 'Ecological theory helps reveal human decision rules from empirical fishery data’

This respiratory contains the simulation and analysis code for the following manuscript:
> Ecological theory helps reveal human decision rules from empirical fishery data
> Yun Ho*, Ke-Yang Chang, Ting-Chun Kuo and Hsi-Cheng Ho

## File structure
The whole analysis is divided into three parts: abundance distribution prediction, individual-based model, and model performance evaluation, each in its own folder. 

The contents in each folder include:
- `01_abundance_predict` performs the model training for GLMM and GAM to generate the daily resource abundance used in simulation.
- `02_IBM` contains IBM of optimal foraging fishers under the four different scenarios.
- `03_model_performance` handles the calculation of model performance metric and the comparison between different model scenarios.

