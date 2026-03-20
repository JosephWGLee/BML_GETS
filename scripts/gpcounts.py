import pandas as pd
import numpy as np
import random
import gpflow 
from IPython.display import display
import tensorflow as tf 

counts_df = pd.read_csv(r"all_normalised_counts.csv", index_col=[0])
sampletable_df = pd.read_csv(r"all_sampletable.csv", index_col=[0])

#CSV processing

counts_df = counts_df.loc[(counts_df.iloc[:, 1:] > 0).sum(axis=1) > 15]
counts_df = counts_df.round().astype(int) 
genes_name = counts_df.index.tolist()

from GPcounts.RNA_seq_GP import rna_seq_gp

X = sampletable_df[["Time"]]
Y = counts_df

#Two sample test

sparse = True
gp_counts = rna_seq_gp(X,Y.loc[genes_name], M=10, safe_mode=False) 
lik_name = "Negative_binomial" 

results = gp_counts.Two_samples_test(lik_name)
results_df = pd.DataFrame(results)

import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

#Statistical Tests

llr_values = results.loc[:,"log_likelihood_ratio"]

degrees_freedom = 1
p_values = [1 - stats.chi2.cdf(abs(llr), degrees_freedom) for llr in llr_values]

results_df_chisq = 1
results_df["p_value"] = 1 - stats.chi2.cdf(abs(results_df["log_likelihood_ratio"]), results_df_chisq)

p_value = results_df["p_value"].values
rejected, q_values, _, _ = multipletests(p_values, method='fdr_bh')
results_df["q_value"] = q_values
results_df["significant"] = rejected

results_df.to_csv("all_gpcounts_significant_genes.csv", index=True)

import os

#Kernel fitting

counts_array = Y.values
counts_array = np.array(Y.values, dtype=np.float64)

time_points_array = np.array(X.values, dtype=np.float64)
time_points_array = time_points_array.reshape(-1,1)

rbf_kernel = gpflow.kernels.SquaredExponential()
periodic_kernel = gpflow.kernels.Periodic(base_kernel=rbf_kernel)
combined_kernel = rbf_kernel + periodic_kernel

#Included due to GPcounts' tendency to fail 

file_path = "Only_RBF_Periodicity_Values.csv"
if os.path.exists(file_path):
    os.remove(file_path)

results = []

for gene_index in range(counts_array.shape[0]):
    
    y = counts_array[gene_index, :].reshape(-1, 1)

    y = tf.convert_to_tensor(y, dtype=tf.float64)

    rbf_value = rbf_kernel(time_points_array, time_points_array)
    periodic_value = periodic_kernel(time_points_array, time_points_array)
    combined_value = rbf_value + periodic_value

    rbf_kernel= gpflow.kernels.SquaredExponential()
    periodic_kernel = gpflow.kernels.Periodic(base_kernel=gpflow.kernels.SquaredExponential())
    combined_kernel = rbf_kernel + periodic_kernel 

    model = gpflow.models.GPR(data=(time_points_array, y), kernel=combined_kernel)
    gpflow.optimizers.Scipy().minimize(model.training_loss, model.trainable_variables)
    
    ratio = rbf_value/combined_value

    results.append({
        "Gene Name" : gene_index,
        "RBF Value": np.mean(rbf_value.numpy()),
        "Peridodicity Value": np.mean(periodic_value.numpy()),
        "Combined Value": np.mean(combined_value.numpy()),
        "Ratio RBF:Periodicity": np.mean(ratio.numpy())
    })

kernel_df = pd.DataFrame(results)
kernel_df.to_csv("all_rbf_periodicity_values.csv", index=False, mode="w")
existing_df = pd.read_csv("all_gpcounts_significant_genes.csv")

updated_df = pd.concat([existing_df, kernel_df.iloc[:, 1:]], axis=1)
updated_df = updated_df.round(3)
updated_df.to_csv("final_gpcounts.csv", index=False)
