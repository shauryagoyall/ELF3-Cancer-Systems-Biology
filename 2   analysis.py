import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import umap
import scipy.stats as ss

"""This file contains the functions that generates a clustermap for the simulated data"""

def plot_heatmap(df,network_name,path_to_plots, replicate):
	
	array_max_min = np.array([df.max().max(),-1*df.min().min()]).max()
	sns.set(font_scale=1.4)
	sns.clustermap(df,cmap='seismic',center = 0, vmin=-1*array_max_min,vmax=array_max_min,annot_kws={"size": 16})
	#plt.title("Heatmap")

	plt.savefig(path_to_plots+network_name+"_heatmap.png",dpi=500)
	plt.close()
	print(replicate + " Heatmap done")

def plot_correlation(df,network_name,path_to_plots, replicate):
	
	sns.set(font_scale=1.4)
	sns.clustermap(df.corr(),cmap='seismic',vmin=-1, center =0, vmax=1, annot_kws={"size": 16})
	#plt.title("Correlation")

	plt.savefig(path_to_plots+network_name+"_correlation.png",dpi=500)
	plt.close()
	print(replicate + " Correlation done")

def multiple_gaussian_fitting(array_to_be_plotted,gene_name,path_to_plots,num_gaussians,sample_size,max_clusts,phenotype):

	AIC_all = {}
	BIC_all = {}

	file = open(path_to_plots+'_'+phenotype+'score_stats.txt','w')
	for i in range(1,max_clusts+1):
		AIC_all[i] = []
		BIC_all[i] = []
	for sample in range(1,sample_size):

		random_state = np.random.RandomState(seed=1)

		X = np.array(array_to_be_plotted).reshape(-1, 1)

		# fit models with 1-10 components
		N = np.arange(1, max_clusts+1)
		models = [None for i in range(len(N))]

		for i in range(len(N)):
			models[i] = GaussianMixture(N[i]).fit(X)

		# compute the AIC and the BIC
		AIC = [m.aic(X) for m in models]
		BIC = [m.bic(X) for m in models]

		for idx,val in enumerate(AIC):
			AIC_all[idx+1].append(AIC[idx])
			BIC_all[idx+1].append(BIC[idx]) 

		M_best = models[num_gaussians-1]              # change number of gaussians here
		s = ""
		mean_array = []
		if num_gaussians>2:
			score_array = np.empty((15,9))
		else:
			score_array = np.empty((15,9))
			
		# for x in range(num_gaussians):
		# 	s += str(M_best.means_[x][0])+"\t"+str(M_best.covariances_[x][0][0])+"\t"+str(M_best.weights_[x])+"\t"
		# 	#print(M_best.means_[0][0],"\t",M_best.means_[1][0],"\t",M_best.means_[2][0],"\t",M_best.covariances_[0][0][0],"\t",M_best.covariances_[1][0][0],"\t",M_best.covariances_[2][0][0],"\t",M_best.weights_[0],"\t",M_best.weights_[1],"\t",M_best.weights_[2])
		# print(s)
		for x in range(num_gaussians):
			mean_array.append(M_best.means_[x][0])

		max_mean = max(mean_array)
		min_mean = min(mean_array)
		y = mean_array.index(max_mean)
		z = mean_array.index(min_mean)

		s += str(M_best.means_[y][0])+"\t"+str(M_best.covariances_[y][0][0])+"\t"+str(M_best.weights_[y])+"\t"
		s += str(M_best.means_[z][0])+"\t"+str(M_best.covariances_[z][0][0])+"\t"+str(M_best.weights_[z])+"\t"
		
		if(num_gaussians>2):
			for k in range(num_gaussians):
				if (k != y and k != z):
					s += str(M_best.means_[k][0])+"\t"+str(M_best.covariances_[k][0][0])+"\t"+str(M_best.weights_[k])+"\t"
					hybrid = M_best.means_[k][0]
		file.write(s+"\n")
		#print(s)

	file.close()

	# print(AIC_all)
	# print(BIC_all)

	fig = plt.figure()
	#fig.subplots_adjust(left=0.12, right=0.97,bottom=0.21, top=0.9, wspace=0.5)

	# plot 1: data + best-fit Mixture
	fig, ax = plt.subplots(figsize=(10,5))

	x = np.linspace(-6, 6, 1000)
	logprob = M_best.score_samples(x.reshape(-1, 1))
	responsibilities = M_best.predict_proba(x.reshape(-1, 1))
	pdf = np.exp(logprob)
	pdf_individual = responsibilities * pdf[:, np.newaxis]

	ax.hist(X, 30, density=True, histtype='stepfilled', alpha=0.4)
	ax.plot(x, pdf, '-k')
	ax.plot(x, pdf_individual, '--k')
	#ax.text(0.04, 0.96, "Best-fit Mixture",ha='left', va='top', transform=ax.transAxes,fontsize=18)
	ax.set_xlim(-2.5,2.5)
	ax.set_xlabel('score', fontsize=16)
	ax.set_ylabel('Frequency', fontsize=16)
	if phenotype == 'EMT' or phenotype == "PDL1":
			plt.axvline(x = hybrid - (abs(hybrid - min_mean)/2), color ='r')
			plt.axvline(x = hybrid + (abs(hybrid - max_mean)/2), color ='r')
	else:
		plt.axvline(x = min_mean +(abs(min_mean - max_mean)/2), color ='r')
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.savefig(path_to_plots+gene_name+".png",dpi=500)
	plt.close()

	# removing the first occurance 
	# AIC_all.pop(1)
	# BIC_all.pop(1)

	AIC_dataframe = pd.DataFrame.from_dict(AIC_all,orient='columns')
	BIC_dataframe = pd.DataFrame.from_dict(BIC_all,orient='columns')

	AIC_dataframe.to_csv(path_to_plots+gene_name+"_AIC.dat",sep="\t",header=None)
	BIC_dataframe.to_csv(path_to_plots+gene_name+"_BIC.dat",sep="\t",header=None)

	# # plot 2: AIC plots
	fig, ax = plt.subplots(figsize=(10,5))
	ax.boxplot(AIC_all.values())
	ax.set_xticklabels(AIC_all.keys())
	plt.xlabel('Cluster Count')
	plt.ylabel('AIC metric')
	plt.title('AIC Values by Cluster Count')

	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.savefig(path_to_plots+gene_name+"_AIC.png",dpi=500)
	plt.close()
	
	# # plot 3: BIC plots
	fig, ax = plt.subplots(figsize=(10,5))
	ax.boxplot(BIC_all.values())
	ax.set_xticklabels(BIC_all.keys())
	plt.xlabel('Cluster Count')
	plt.ylabel('BIC metric')
	plt.title('BIC Values by Cluster Count')

	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.savefig(path_to_plots+gene_name+"_BIC.png",dpi=500)
	plt.close()
    
def plot_histograms(df, network_name, plot_path, replicate):
    
	#df["EM_score"] = (df['ZEB'] + df['SLUG'] - df['CDH1'] - df['MIR'])/4
	#multiple_gaussian_fitting(df["EM_score"],"EMT_score_dist",plot_path,3,15,6,'EMT')
	#print(replicate + " EM Score Histogram done")
	
	#df["EM_score"] = (df['ZEB'] + df['SLUG'] - df['CDH1'] - df['MIR'])/4
	#multiple_gaussian_fitting(df["PDL1"],"PDL1_dist",plot_path,3,15,6,'PDL1')
	#print(replicate + " PDL1 Score Histogram done")
	
	df["TamRes_score"] = df["ERA36"] - df["ERA66"]
	multiple_gaussian_fitting(df["TamRes_score"],"TamRes_score_dist",plot_path,2,15,6,'TamRes')
	print(replicate + " TamRes Score Histogram done")
 
#This function performs the Umap analysis. 
def UMAP_analysis(data,n_neighbors=15, min_dist=0.1, n_components=2, metric='euclidean'):
    fit = umap.UMAP(n_neighbors=n_neighbors,min_dist=min_dist,n_components=n_components,metric=metric)
    umap_data = fit.fit_transform(data)

    return umap_data

#This function plots the output from the Umap analysis. 
def plot_UMAP_scatter(df, plot_path, replicate):
    sub_dataframe = df.sample(n=10000)
    embedding = UMAP_analysis(sub_dataframe,n_neighbors=100)
    UMAP = pd.DataFrame()

    UMAP['UMAP_1'] = embedding[:,0]
    UMAP['UMAP_2'] = embedding[:,1]
    
    #### UMAP of EMT score
    
    UMAP.plot.scatter(x = 'UMAP_1',y='UMAP_2',s=1,c=((sub_dataframe['ZEB'] + sub_dataframe['SLUG'] - sub_dataframe['CDH1'] - sub_dataframe['MIR'])/4),cmap=plt.cm.nipy_spectral)
    plt.title('EMT score UMAP',fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(plot_path+"EM_score_Umap.png",dpi=400)
    plt.close()
    print(replicate + " EM score UMAP done")

    #### UMAP of tamoxifen resistance score
    
    # UMAP.plot.scatter(x='UMAP_1',y='UMAP_2',s=1,c=sub_dataframe['ERa36']-sub_dataframe['ERa66'],cmap=plt.cm.nipy_spectral)
    # plt.title('Tamoxifen resistance score UMAP',fontsize=14)
    # plt.xlabel('UMAP_1',fontsize=14)
    # plt.ylabel('UMAP_2',fontsize=14)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.savefig(plot_path+"TamRes_Umap=.png",dpi=400)
    # plt.close()   
 
# This function reads the file containing the mean and standard deviation for each gaussian in the score distribution.
def read_scoring(plot_path,phenotype):

	data = pd.read_csv(plot_path+'_'+phenotype+'score_stats.txt',delimiter="\t",header=None)
	mean_of_the_columns = data.mean(axis = 0)
	sdev_of_the_columns = data.std(axis = 0)

	return mean_of_the_columns, sdev_of_the_columns

# This function uses mean and standard deviation for each gauassian in the score distribution determine a boundary condition
def phenotypes_boundary(mean_EMT,std_EMT):
	boundary = []
	boundary.append(mean_EMT[6] -(abs(mean_EMT[6]-mean_EMT[3])/2))
	boundary.append(mean_EMT[6] +(abs(mean_EMT[6]-mean_EMT[0])/2))
	return boundary

# This function generates scatter plot of EM score vs PD-L1 expression
def plot_EM_PD_scatter(df, plot_path, replicate):

	mean_EMT, std_EMT = read_scoring(plot_path,'EMT')
	boundaries = phenotypes_boundary(mean_EMT,std_EMT) 
    
	df['EM_score'] = (df['ZEB'] + df['SLUG'] - df['CDH1'] - df['MIR'])/4
	plt.scatter(df['EM_score'], df['PDL1'],marker="o",s=0.1,c='black')
	plt.axvline(x=boundaries[0], color='r')
	plt.axvline(x=boundaries[1], color='r')
	plt.axhline(y=0, color='r')
	plt.xlabel('EMT_score')
	plt.ylabel('PD-L1 Expression')
	plt.tight_layout()
	plt.savefig(plot_path + "EM_PD_scatter.png" , dpi=500)
	plt.close()
	print(replicate + " EM-PD Scatter done")
    
	print(ss.spearmanr(df["EM_score"],df["PDL1"]))
    
# This function generates scatter plot of EM score vs PD-L1 expression
def plot_TamRes_PD_scatter(df, plot_path, replicate):

	df['TamRes_score'] = df["ERA36"] - df["ERA66"]
	plt.scatter(df['TamRes_score'], df['PDL1'],marker="o",s=0.1,c='black')
	plt.axvline(x=0, color='r')
	plt.axhline(y=0, color='r')
	plt.xlabel('Tamoxifen Resistance')
	plt.ylabel('PD-L1 Expression')
	plt.tight_layout()
	plt.savefig(plot_path + "TamRes_PD_scatter.png" , dpi=500)
	plt.close()
	print(replicate + " TamRes-PD Scatter done")
    
	print(ss.spearmanr(df["TamRes_score"],df["PDL1"]))
    
def plot_TamRes_PD_tri_scatter(df, plot_path, replicate):

	mean_EMT, std_EMT = read_scoring(plot_path,'PDL1')
	boundaries = phenotypes_boundary(mean_EMT,std_EMT) 
    
	df['TamRes_score'] = df["ERA36"] - df["ERA66"]
	plt.scatter(df['TamRes_score'], df['PDL1'],marker="o",s=0.1,c='black')
	plt.axhline(y=boundaries[0], color='r')
	plt.axhline(y=boundaries[1], color='r')
	plt.axvline(x=0, color='r')
	plt.xlabel('Tamoxifen Resistance')
	plt.ylabel('PD-L1 Expression')
	plt.tight_layout()
	plt.savefig(plot_path + "TamRes_PD_tri_scatter.png" , dpi=500)
	plt.close()
	print(replicate + " EM-PD Tri Scatter done")
    
# This function generates plot of conditional probability of PD-L1 high given a phenotype
def plot_prob_state(df, plot_path, replicate):

	mean_EMT, std_EMT = read_scoring(plot_path,'EMT')
	boundaries = phenotypes_boundary(mean_EMT,std_EMT) 
    
	df['EM_score'] = (df['ZEB'] + df['SLUG'] - df['CDH1'] - df['MIR'])/4
    ## state = [pdl1 low, pdl1 high, prob]
	E = [0,0,0]
	H = [0,0,0]
	M = [0,0,0]
	for i in range(len(df['EM_score'])):
	    val = df['EM_score'][i]
	    if val <= boundaries[0]:
	        if df['PDL1'][i] <= 0:
	            E[0] +=1
	        else:
	            E[1] +=1
	    elif val > boundaries[0] and val <= boundaries[1]:
	        if df['PDL1'][i] <= 0:
	            H[0] +=1
	        else:
	            H[1] +=1
	    else:
	        if df['PDL1'][i] <= 0:
	            M[0] +=1
	        else:
	            M[1] += 1
                
	E[2]= E[1] / ( E[0] + E[1] )
	H[2]= H[1] / ( H[0] + H[1] )
	M[2]= M[1] / ( M[0] + M[1] )
    
	prob = [E[2], H[2], M[2]]
    
	
	x = np.arange(3)
	plt.bar(x, prob, color = ['black', 'orange','red' ], alpha=0.6)
	plt.ylim(0, 1) 
	plt.xticks(x, ['Epithelial','Hybrid','Mesenchymal'] )
	plt.xlabel("Phenotype", weight = 'bold')
	plt.ylabel("Conditional Probability", weight = 'bold')
	plt.tight_layout()
	plt.savefig(plot_path + "prob_state.png" , dpi=500)
	plt.close()
	print(replicate + " prob vs phenotype done")

# This function generates plot of conditional probability of PD-L1 high given a phenotype with PDL1 3 levels
def plot_prob_state_tri(df, plot_path, replicate):

	mean_EMT, std_EMT = read_scoring(plot_path,'EMT')
	boundaries = phenotypes_boundary(mean_EMT,std_EMT)
    
	mean_PD, std_PD = read_scoring(plot_path,'PDL1')
	boundaries_PD = phenotypes_boundary(mean_PD,std_PD) 
    
	df['EM_score'] = (df['ZEB'] + df['SLUG'] - df['CDH1'] - df['MIR'])/4
    ## state = [pdl1 low, pdl1 med, pdl1 high, prob]
	E = [0,0,0,0]
	H = [0,0,0,0]
	M = [0,0,0,0]
	for i in range(len(df['EM_score'])):
	    val = df['EM_score'][i]
	    if val <= boundaries[0]:
	        if df['PDL1'][i] <= boundaries_PD[0]:
	            E[0] +=1
	        elif df['PDL1'][i] > boundaries_PD[0] and df['PDL1'][i] <= boundaries_PD[1]:
	            E[1] +=1
	        else :
	            E[2] +=1
	    elif val > boundaries[0] and val <= boundaries[1]:
	        if df['PDL1'][i] <= boundaries_PD[0]:
	            H[0] +=1
	        elif df['PDL1'][i] > boundaries_PD[0] and df['PDL1'][i] <= boundaries_PD[1]:
	            H[1] +=1
	        else :
	            H[2] +=1
	    else:
	        if df['PDL1'][i] <= boundaries_PD[0]:
	            M[0] +=1
	        elif df['PDL1'][i] > boundaries_PD[0] and df['PDL1'][i] <= boundaries_PD[1]:
	            M[1] +=1
	        else :
	            M[2] +=1
                
	E[3]= E[2] / ( E[0] + E[1] + E[2] )
	H[3]= H[2] / ( H[0] + H[1] + H[2] )
	M[3]= M[2] / ( M[0] + M[1] + M[2] )
    
	prob = [E[3], H[3], M[3]]
    
	
	x = np.arange(3)
	plt.bar(x, prob, color = ['black', 'orange','red' ], alpha=0.6)
	plt.ylim(0, 1) 
	plt.xticks(x, ['Epithelial','Hybrid','Mesenchymal'] )
	plt.xlabel("Phenotype", weight = 'bold')
	plt.ylabel("Conditional Probability", weight = 'bold')
	plt.tight_layout()
	plt.savefig(plot_path + "prob_state_tri.png" , dpi=500)
	plt.close()
	print(replicate + " prob vs phenotype tri done")
    
# This function generates plot of conditional probability of PD-L1 high given a resistance
def plot_prob_resistance(df, plot_path, replicate):

	df['TamRes_score'] = df["ERA36"] - df["ERA66"]
    ## resistance level = [pdl1 low, pdl1 high, prob]
	S = [0,0,0]
	R = [0,0,0]

	for i in range(len(df['TamRes_score'])):
	    val = df['TamRes_score'][i]
	    if val <= 0:
	        if df['PDL1'][i] <= 0:
	            S[0] +=1
	        else:
	            S[1] +=1
	    else:
	        if df['PDL1'][i] <= 0:
	            R[0] +=1
	        else:
	            R[1] += 1
                
	S[2]= S[1] / ( S[0] + S[1] )
	R[2]= R[1] / ( R[0] + R[1] )
    
	prob = [S[2], R[2]]
    
	x = np.arange(2)
	plt.bar(x, prob, color = ['green', 'maroon'], alpha=0.6)
	plt.ylim(0, 1) 
	plt.xticks(x, ['Sensitive','Resistant'] )
	plt.xlabel("Phenotype", weight = 'bold')
	plt.ylabel("Conditional Probability", weight = 'bold')
	plt.tight_layout()
	plt.savefig(plot_path + "prob_resistance.png" , dpi=500)
	plt.close()
	print(replicate + " prob vs resistance done")

# This function generates plot of conditional probability of PD-L1 high given a resistance
def plot_prob_resistance_tri(df, plot_path, replicate):

	df['TamRes_score'] = df["ERA36"] - df["ERA66"]
    ## resistance level = [pdl1 low, pdl1 med, pdl1 high, prob]
	S = [0,0,0,0]
	R = [0,0,0,0]
    
	mean_PD, std_PD = read_scoring(plot_path,'PDL1')
	boundaries_PD = phenotypes_boundary(mean_PD,std_PD)
    
	for i in range(len(df['TamRes_score'])):
	    val = df['TamRes_score'][i]
	    if val <= 0:
	        if df['PDL1'][i] <= boundaries_PD[0]:
	            S[0] +=1
	        elif df['PDL1'][i] > boundaries_PD[0] and df['PDL1'][i] <= boundaries_PD[1]:
	            S[1] +=1
	        else :
	            S[2] +=1
	    else:
	        if df['PDL1'][i] <= boundaries_PD[0]:
	            R[0] +=1
	        elif df['PDL1'][i] > boundaries_PD[0] and df['PDL1'][i] <= boundaries_PD[1]:
	            R[1] +=1
	        else :
	            R[2] +=1
                
	S[3]= S[2] / ( S[0] + S[1] +S[2] )
	R[3]= R[2] / ( R[0] + R[1] +R[2])
    
	prob = [S[3], R[3]]
    
	x = np.arange(2)
	plt.bar(x, prob, color = ['green', 'maroon'], alpha=0.6)
	plt.ylim(0, 1) 
	plt.xticks(x, ['Sensitive','Resistant'] )
	plt.xlabel("Phenotype", weight = 'bold')
	plt.ylabel("Conditional Probability", weight = 'bold')
	plt.tight_layout()
	plt.savefig(plot_path + "prob_resistance_tri.png" , dpi=500)
	plt.close()
	print(replicate + " prob vs resistance tri done")
    
    
replicates = ['r1','r2']

def all_analysis(replicate):
	core_path = "./elf3/"
	
	network_name = "elf3"
	data_path = core_path + replicate + "/"
	#plot_path = data_path + "plots/"
    ## For 3 levels of PDL1 use below path 
	plot_path = data_path + "tri_plots/"
	if not os.path.exists(plot_path):
         os.makedirs(plot_path)
	df = pd.read_excel(data_path + network_name +".xlsx")
    
    ##############################################
    ##############################################
    ##############      Plots       ##############
    ##############################################
    ##############################################
    
	#plot_heatmap(df, network_name, plot_path , replicate)
	#plot_correlation(df,network_name,plot_path, replicate)
	#plot_histograms(df, network_name, plot_path, replicate)
	#plot_EM_PD_scatter(df, plot_path,  replicate)
	#plot_TamRes_PD_scatter(df, plot_path,  replicate)
	#plot_prob_state(df, plot_path, replicate)
	#plot_prob_resistance(df, plot_path, replicate)
	#plot_UMAP_scatter(df, plot_path, replicate)
	#plot_prob_state_tri(df, plot_path, replicate)
	#plot_TamRes_PD_tri_scatter(df, plot_path, replicate)
	plot_prob_resistance_tri(df, plot_path, replicate)
	
for replicate in replicates:
	all_analysis(replicate)