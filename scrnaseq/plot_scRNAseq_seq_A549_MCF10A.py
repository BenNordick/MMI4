import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

filea = "A549.csv" #GSE147405 Cook and Vanderhyden 2020
filea = "MCF10A.csv" #GSE213753 Panchy et al. 2022
dfa = pd.read_csv(filea, index_col=0)
dfa['200b_d'] = dfa.TS_200b - dfa.background
dfa['101p1_d'] = dfa.TS_101p1 - dfa.background
dfa['CDH1_d'] = dfa.CDH1 - dfa.background
dfa['SNAI1_d'] = dfa.SNAI1 - dfa.background
dfa['SNAI2_d'] = dfa.SNAI2 - dfa.background
dfa['ZEB1_d'] = dfa.ZEB1 - dfa.background
dfa['VIM_d'] = dfa.VIM - dfa.background
print(dfa.columns)

if filea.startswith("A549"):
    hue_order = ["0d", "8h", "1d", "3d", "7d"]#, "8h_rm", '1d_rm', "3d_rm"]
elif filea.startswith("MCF10A"):
    hue_order = [0, 12.5, 25, 50, 100, 200]
cset = sns.color_palette("husl", 8)
cset2 = sns.color_palette("Spectral", 6)
palette = {0:cset[0], 12.5: cset[1], 25: cset[2], 
        50:cset[3], 100:cset[4], 200:cset[5],
        "0d":cset[0], "8h": cset[1], "1d": cset[2], 
        "3d":cset[3], "7d":cset[4]
        }

if filea.startswith("A549"):
    hue = "Time"
    xv, xl = "Time", "Time (days)"
elif filea.startswith("MCF10A"):
    xv, xl = "GBC_pM", "Dose (pM)"

fig, ax = plt.subplots(figsize=(5,4))
fig.subplots_adjust(left=0.15, bottom=0.135, right=0.75)
ax = sns.scatterplot(x="UMAP_1", y="UMAP_2", data=dfa, hue=hue, ax=ax, alpha=0.3, size=1, hue_order=hue_order, palette=palette)#, legend=False, palette=palette)
ax.legend(bbox_to_anchor=(1.1, 0, 0.3, 0.9))
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
plt.show()

fig, ax = plt.subplots(figsize=(4,4))
fig.subplots_adjust(left=0.25, bottom=0.135)
sns.boxplot(x=xv, y="101p1_d", data=dfa, order=hue_order, palette=palette)#, hue="Time", ax=ax, alpha=0.3, size=1, hue_order=hue_order)#, legend=False, palette=palette)
#sns.boxplot(x=xv, y="200b_d", data=dfa, order=hue_order, palette=palette)#, hue="Time", ax=ax, alpha=0.3, size=1, hue_order=hue_order)#, legend=False, palette=palette)
#sns.boxplot(x=xv, y="CDH1_d", data=dfa, order=hue_order, palette=palette)#, hue="Time", ax=ax, alpha=0.3, size=1, hue_order=hue_order)#, legend=False, palette=palette)
#sns.boxplot(x=xv, y="CDH1_d", data=dfa[dfa.CDH1_d>0.05], order=hue_order, palette=palette)#, hue="Time", ax=ax, alpha=0.3, size=1, hue_order=hue_order)#, legend=False, palette=palette)
#sns.boxplot(x=xv, y="ZEB1_d", data=dfa[dfa.ZEB1_d>0.5], order=hue_order, palette=palette)#, hue="Time", ax=ax, alpha=0.3, size=1, hue_order=hue_order)#, legend=False, palette=palette)
#sns.boxplot(x=xv, y="VIM_d", data=dfa, order=hue_order, palette=palette)#, hue="Time", ax=ax, alpha=0.3, size=1, hue_order=hue_order)#, legend=False, palette=palette)
#sns.boxplot(x=xv, y="SNAI1_d", data=dfa[dfa.SNAI1_d>0.5], order=hue_order, palette=palette)#, hue="Time", ax=ax, alpha=0.3, size=1, hue_order=hue_order)#, legend=False, palette=palette)
#sns.boxplot(x=xv, y="SNAI2_d", data=dfa[dfa.SNAI2_d>-100*0.5], order=hue_order, palette=palette)#, hue="Time", ax=ax, alpha=0.3, size=1, hue_order=hue_order)#, legend=False, palette=palette)
ax.set_ylabel('Mean miR-101 Target\nGene Expression (Z-score)')
#ax.set_ylabel('Mean miR-200bc Target\nGene Expression (Z-score)')
#ax.set_ylabel('CDH1 Expression (Z-score)')
#ax.set_ylabel('CDH1 Expression* (Z-score)')
#ax.set_ylabel('ZEB1 Expression* (Z-score)')
#ax.set_ylabel('VIM Expression (Z-score)')
#ax.set_ylabel('SNAI1 Expression* (Z-score)')
#ax.set_ylabel('SNAI2 Expression (Z-score)')
ax.set_xlabel(xl)
plt.show()

