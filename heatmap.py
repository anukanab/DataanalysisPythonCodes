import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df=pd.read_csv('/data/salomonis2/NCI-R01/hg19-signatures/concordance.txt', sep='\t', index_col=0)

#heatmap=sns.heatmap(df, annot=True, yticklabels=False)
sns.heatmap(df, annot=True)
plt.savefig("/data/salomonis2/NCI-R01/hg19-signatures/concordance_heatmap.png")

#sns.plt.show()