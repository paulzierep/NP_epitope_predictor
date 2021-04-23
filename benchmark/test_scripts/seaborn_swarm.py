import seaborn as sns
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1,1,figsize = (11,6))

sns.set(style="whitegrid")
tips = sns.load_dataset("tips")

print(tips)

sns.swarmplot(x=tips["total_bill"], ax = ax)

plt.show()