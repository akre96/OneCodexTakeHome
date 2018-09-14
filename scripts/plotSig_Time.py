import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

NUM_TO_ANALYZE=3

# LOAD DATA
SigGenus_List= pd.read_csv('../results/SignificantGenus_P05_TimeChange.csv').sort_values('P-Value')
SigGenus_Reads= pd.read_csv('../results/SignificantGenus_P05_TimeChange_Reads.csv')

# Create list of genuses to plot
GenusToAnalyze=SigGenus_List['Genus'][:NUM_TO_ANALYZE].values.tolist()
GenusToAnalyze.append('Bifidobacterium')

plotDF=SigGenus_Reads.loc[SigGenus_Reads['Genus'].isin(GenusToAnalyze)]


# Plot

ax = sns.lineplot(x='time',y='Read Counts',hue='Genus',size='P-Value',data=plotDF)
ax.set(xlabel="Time (Days)")

plt.title("TAC Mice Galvanized with JAX")
plt.show()



