######################################################################################
#                   Generate mean |SHAP| for each architecture

# This is a script for taking mean of SHAP values and concatinating results for each
# fold in every architecture, preparing the result for plotting. The script assumes that
# SHAP-data has been prepared for each fold using 'sim_shapley.py' script used on Markov.
# 'dset = "70k"' default, but "180k" possible.
######################################################################################

#Libraries
import numpy as np
import pandas as pd
import pyarrow.feather as feather

# Set architecture to run
arch = 2
# Choose 70k or 180k dataset
dset = "70k"

shap_list = []
for i in range(1,11):
    print("Starting run ", str(i))
    name = f'C:/Users/gard_/Documents/MasterThesis/Code/Results/{dset}/SHAP/shap_{dset}_arch_{arch}_fold_{str(i)}.csv'
    temp_shap = pd.read_csv(name)
    shap_list.append(temp_shap)

mean_shap = np.abs(np.concatenate(shap_list, axis=0)).mean(axis=0) # (p x 1), mean over the indvs

if mean_shap.ndim == 2:
    mean_shap = mean_shap.flatten()

# Save to be further manipulated in Chromosome2Manhattan.R
shap_df = pd.DataFrame(mean_shap.reshape(1,-1), columns=temp_shap.columns) # (p x 1), Values for Manhattan

feather.write_feather(shap_df, f'C:/Users/gard_/Documents/MasterThesis/Code/Results/{dset}/SHAP/shap_{dset}_arch_{arch}.feather')

print(f'Done with architecture {arch}!')