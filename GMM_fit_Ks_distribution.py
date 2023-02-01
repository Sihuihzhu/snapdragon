import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pandas import DataFrame
from sklearn import datasets
import scipy.stats as ss
from sklearn.mixture import GaussianMixture
import seaborn as sns
import sys
import matplotlib.font_manager as font_manager
font = font_manager.FontProperties(style='italic', size=14)

def get_array_for_mixture(df):
	X = np.array(df["Ks"].dropna())
	X = X[X > 0]
	X = np.log(X).reshape(-1, 1)
	return X

def plot_mixture(model, data, ax, l=0, u=5, color='black', label='',alpha=0.2,log=False, bins=25, alpha_l1=1):
	x = np.linspace(l, u, 1000).reshape((-1, 1))
	if not log:
		data = np.exp(data)
	data = data[data >= l]
	data = data[data <= u]
	means = model.means_
	varcs = model.covariances_
	weights = model.weights_
	mix = None
	first = True
	for k in range(len(means)):
		if not log:
			curve = ss.lognorm.pdf(x, scale=np.exp(means[k]), s=np.sqrt(varcs[k])) * weights[k]
		else:
			curve = ss.norm.pdf(x, loc=means[k], scale=np.sqrt(varcs[k])) * weights[k]
		if first:
			mix = curve
			first = False
		else:
			mix += curve
	ax.plot(x, mix, color=color,label=label, alpha=alpha_l1)
	ax.set_xlim(l, u)
	ax.set_xlabel("$K_{\mathrm{S}}$")
	ax.set_ylabel("Density")
	return ax

def filter_group_data(df, aln_id=0, aln_len=300, aln_cov=0, min_ks=0.5, max_ks=5,weights_outliers_included=False):
	df = df.dropna()
	df = df[df["AlignmentCoverage"] >= aln_cov]
	df = df[df["AlignmentIdentity"] >= aln_id]
	df = df[df["AlignmentLength"] >= aln_len]
	if not weights_outliers_included:
		df = df[df["Ks"] > min_ks]
		df = df[df["Ks"] <= max_ks]
	df = df.groupby(['Family', 'Node']).mean()
	if weights_outliers_included:
		df = df[df["Ks"] > min_ks]
		df = df[df["Ks"] <= max_ks]
	return df


plot_order =['Antirrhinum_hispanicum','Misopates_orontium','Collinsia_rattanii','Salvia_miltiorrhiza']

component_d= {
				'Antirrhinum_hispanicum':[2,1.5],
				'Collinsia_rattanii':[2,1.5],
				'Salvia_miltiorrhiza':[2,1.5],
				'Misopates_orontium':[2,1.0]}
mcl_file = {'Antirrhinum_hispanicum':'Antirrhinum_hispanicum_ksd/Antirrhinum_hispanicum.cds.ks.tsv','Collinsia_rattanii':'wgd_ksd/Collinsia_rattanii.cds.ks.tsv','Salvia_miltiorrhiza':'Salvia_miltiorrhiza_ksd/Salvia_miltiorrhiza.cds.ks.tsv','Misopates_orontium':'Misopates_orontium_ksd/Misopates_orontium.cds.ks.tsv'}
species = plot_order 
fig, ax = plt.subplots(1,1 , figsize=(6, 6))
colors = ['C1','C2','C3','C4']
ks_filter = [0.2,0.2,0.3,0.3]
for i in range(len(species)):
	print species[i]
	df = pd.read_csv(mcl_file[species[i]],index_col=0, sep='\t')
	df = filter_group_data(df, 0,300,0.5,ks_filter[i], 3)
	X = get_array_for_mixture(df)
	if i==2:X=X-0.1
	models=[None]
	models[0] = GaussianMixture(n_components=component_d[species[i]][0],covariance_type='full', max_iter=100,n_init=1).fit(X)
	data = X 
	plot_mixture(models[0], data, ax, l=0, u=3, bins=100,color=colors[i],label=species[i].replace("_"," "))
	ax.set_ylim(0,component_d[species[i]][1])
	if i==0:
		ax.set_ylabel("Density",fontsize=20)
	else:
		ax.set_ylabel("",fontsize=2)
	ax.set_xlabel("$K_{\mathrm{S}}$",fontsize=20)
ax.set_ylim(0,1.1)
ax.set_xticks([0,1,2,3])
ax.set_xticklabels(['0','1','2','3'])
ax.set_ylabel("Density",fontsize=20)
ax.legend(fancybox=True,frameon=True,prop=font)
sns.despine(offset=5)
fig.tight_layout()
plt.savefig('Ks.png',dpi=600,format="png")
plt.savefig('Ks.svg',dpi=600,format="svg")


