import scipy.stats as stats
import numpy as np
import pandas as pd

p = [0,0,1,1,0,0,0,0,0,0,]
n = [0,0,1,1,0,0,1,1,2,2,2,10,10,10,3,3,3,3,3,3,3,3,3]

def fisher_exact_for_coninuous_outcome(p, n):

	p = np.array(p)
	n = np.array(n)

	p_has_feature = np.count_nonzero(p)
	p_no_feature = np.count_nonzero(p == 0)

	n_has_feature = np.count_nonzero(n)
	n_no_feature = np.count_nonzero(n == 0)

	oddsratio, pvalue = stats.fisher_exact([[p_has_feature, n_has_feature], [p_no_feature, n_no_feature]])

	return(pvalue)

# p_val = fisher_exact_for_coninuous_outcome(p,n)
# print(p_val)

def chi2_undefined_categorical_observation(p, n):

	p = pd.Series(p, name = 'observed')
	n = pd.Series(n, name = 'expected')

	p_counts = p.value_counts()
	n_counts = n.value_counts()

	# print(p_counts)
	# print(n_counts)

	final_df = pd.concat([p_counts, n_counts], axis = 1).T

	final_df = final_df.fillna(0)

	# print(final_df)

	chi2, p_val, dof, expected = stats.chi2_contingency(final_df)

	return(p_val)

p_val = chi2_undefined_categorical_observation(p,n)
# print(p_val)
# oddsratio, p = stats.fisher_exact(final_df)
