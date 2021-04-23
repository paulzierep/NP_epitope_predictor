import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np

plt.rcParams.update({'font.size': "14"}) #bigger text fro phd

###########################
#Dataset choice
###########################

data_storage = os.path.join('benchmark_data')

data_to_use = {}
for folder in os.listdir(data_storage):
    if len(folder) == 11 and not 'MR' in folder:
        
        for x in range(10):
            if str(x) in folder:
                radius = x

        if 'CF' in folder:
            chiral = 'NC FP'
        if 'CT' in folder:
            chiral = 'C FP'

        data_to_use[folder] = '{0}, R.: {1}'.format(chiral, radius)

print(data_to_use)

cell_types = {'b_cell':'B cell','t_cell':'T cell'}

#get some colors for the plot
colors = {}
c_choice = cm.rainbow(np.linspace(0, 1, len(data_to_use)))
for index, data_name in enumerate(data_to_use):
    colors[data_name] = c_choice[index]

###########################
#Plot individual
###########################

# fig1, axs1 = plt.subplots(4,2,figsize = (15,15))
# fig2, axs2 = plt.subplots(4,2,figsize = (15,15))
# axis1 = axs1.reshape(-1)
# axis2 = axs2.reshape(-1)

# base_range = np.around(np.geomspace(1, 2048, num=12)).astype(int) #the reange for all clf even though they might not go that far

# #locate all benchmark files and store cluster/cell type
# for folder in os.listdir(data_storage):
#     if folder in data_to_use:
#         folder_path = os.path.join(data_storage, folder)
#         for sub_folder in os.listdir(folder_path):
#             if sub_folder in cell_types: 
#                 cell_type = sub_folder
#                 sub_folder_path = os.path.join(data_storage, folder, sub_folder)
#                 for benchmark_file in os.listdir(sub_folder_path):
#                     cluster = int(benchmark_file.replace('.csv',''))

#                     print(folder)
#                     print(sub_folder)
#                     print(benchmark_file)
#                     print(cell_type)
#                     print(cluster)
#                     benchmark_file_path = os.path.join(data_storage, folder, sub_folder, benchmark_file)
#                     df = pd.read_csv(benchmark_file_path, index_col=0)

#                     # print(df['mean'])
#                     # exit()

#                     label = "{0}, {1}".format(data_to_use[folder], cell_types[cell_type])
#                     #plot

#                     if cell_type == 'b_cell':
#                         # color = 'blue'
#                         axis = axis1
#                     if cell_type == 't_cell':
#                         # color = "orange"
#                         axis = axis2

#                     color = colors[folder]

#                     # print(df['mean'])
#                     # print(df['err'])
#                     # print(color)

#                     df['mean'].plot(yerr = df['err'], ax = axis[cluster], label=label)

#                     #axis properties    
#                     axis[cluster].set_xscale('log', basex=2)
#                     axis[cluster].set_xticks(base_range)
#                     axis[cluster].set_xticklabels(base_range)
                    
#                     axis[cluster].set_xlabel("Features (chi2)")
#                     axis[cluster].set_ylabel("ROC-AUC") 
                    
#                     # axis[index].axhline(y=sim_data_storage_b[index], ls = "--", color = "blue", label = 'Sim (B cell)')
                    
#                     axis[cluster].axhline(y=0.8, ls = "-", color = "black")
                    
#                     axis[cluster].set_title("Cluster {0}".format(cluster))
#                     axis[cluster].legend(loc = 'lower right')

# fig1.tight_layout()
# fig2.tight_layout()
# fig1.savefig('Chiral_FP_and_radius_vs_AUC_b_cell.png', bbox_inches = "tight")
# fig2.savefig('Chiral_FP_and_radius_vs_AUC_t_cell.png', bbox_inches = "tight")
# exit()


#############################
#Get max AUC score for each query
#############################

main_df = pd.DataFrame()
# main_df.set_index(['cell_type','cluster','data'], inplace=True)

#locate all benchmark files and store cluster/cell type
idx = 0
for folder in os.listdir(data_storage):
    if folder in data_to_use:
        folder_path = os.path.join(data_storage, folder)
        for sub_folder in os.listdir(folder_path):
            if sub_folder in cell_types: 
                cell_type = sub_folder
                sub_folder_path = os.path.join(data_storage, folder, sub_folder)
                for benchmark_file in os.listdir(sub_folder_path):
                    cluster = int(benchmark_file.replace('.csv',''))

                    print(folder)
                    print(sub_folder)
                    print(benchmark_file)
                    print(cell_type)
                    print(cluster)
                    benchmark_file_path = os.path.join(data_storage, folder, sub_folder, benchmark_file)
                    df = pd.read_csv(benchmark_file_path, index_col=0)


                    # print()
                    # max_mean = df.loc[8,'mean']

                    max_mean = df['mean'].max()
                    max_err = df.loc[df['mean'].idxmax(),'err']

                    main_df.loc[idx, 'max_mean'] = max_mean
                    main_df.loc[idx, 'max_err'] = max_err
                    main_df.loc[idx, 'cell_type'] = cell_types[cell_type]
                    main_df.loc[idx, 'cluster'] = cluster
                    main_df.loc[idx, 'data'] = data_to_use[folder] 

                    idx += 1




# print(mai.loc[:,['max_mean','max_err']].mean())
# exit()

#mean for all cluster and cell_type

fig, ax = plt.subplots(1,1,figsize = (8,5))

print(main_df)

mean_per_data = main_df.groupby('data').apply(lambda group: group.loc[:,['max_mean','max_err']].mean())

mean_per_data = mean_per_data.reset_index()

new = mean_per_data['data'].str.split(",", n = 1, expand = True) 

mean_per_data['Chiral'] = new[0]
mean_per_data['Radius'] = new[1]

mean_per_data['Radius'] = mean_per_data['Radius'].str.strip('R.: ')

mean_per_data.to_csv('Chiral_FP_and_radius_vs_AUC.csv')
exit()

mean_per_data.plot.bar(ax = ax, y='max_mean', rot=90, yerr = mean_per_data['max_err'], color=c_choice)

# ind = np.arange(0,mean_per_data.shape[0]/10,0.1)  # the x locations for the groups
# width = 0.1       # the width of the bars

# ax.bar(ind, mean_per_data['max_mean'], width, yerr =  mean_per_data['max_err'])
# ax.set_xticklabels(mean_per_data.index, rotation=45, ha='right')

for p in ax.patches:
    ax.annotate(str(round(p.get_height(),2)), (p.get_x() * 1.005, p.get_height() * 0.90), rotation=90)

ax.get_legend().remove()
ax.set_ylabel('Mean AUC for all \n clusters and epitope categories')
ax.set_xlabel('Morgan fingerprint parameter set')
ax.set_ylim(0.45,1.05)

plt.savefig('Chiral_FP_and_radius_vs_AUC.png', bbox_inches = "tight")
# fig.tight_layout()
# plt.show()

#per cell type
# fig, axs = plt.subplots(1,2,figsize = (15,15))

# for idx, (index, cell_group) in enumerate(main_df.groupby('cell_type')):

# 	mean_per_data = cell_group.groupby('data').apply(lambda group: group['max_mean'].mean())

# 	ax = axs[idx]
# 	mean_per_data.plot.bar(ax = ax)

# 	for p in ax.patches:
# 	    ax.annotate(str(round(p.get_height(),2)), (p.get_x() * 1.005, p.get_height() * 1.005))

# plt.tight_layout()
# plt.show()