from filecmp import dircmp
import os
import numpy as np
import numpy.linalg as la


def relative_entropy(motif, predicted_motif):
    """
    @param motif: path to motif.txt
    @param predicted_motif: path to predictedmotif.txt
    @return entropy: average entropy between motif and predicted_motif
    """
    
    motif_lines = []
    predicted_motif_lines = []
    myMotif = []
    myPredictedMotif = []
    entropy = 0 
    
    with open(motif, 'r') as motif_file:
        motif_lines = motif_file.readlines()
    with open(predicted_motif, 'r') as predicted_motif_file:
        predicted_motif_lines = predicted_motif_file.readlines()
    
    for i in range(1, len(motif_lines)-1):
        myMotif.append(motif_lines[i].strip().split("\t"))
        for j in range(len(myMotif[i-1])):
            myMotif[i-1][j] = float(myMotif[i-1][j].strip())
    
    for i in range(1, len(predicted_motif_lines)-1):
        myPredictedMotif.append(predicted_motif_lines[i].strip().split("\t"))
        for j in range(len(myPredictedMotif[i-1])):
            myPredictedMotif[i-1][j] = float(myPredictedMotif[i-1][j].strip())
    myMotif = np.array(myMotif)
    myPredictedMotif = np.array(myPredictedMotif)
    for real, pred in zip(myMotif, myPredictedMotif):
        real[real==0] = 1e-6
        pred[pred==0] = 1e-6
        for i in range(4): 
            #print(real)
            entropy += pred[i] * (np.log(pred[i]/real[i]))
    
    return entropy


def sites_norm(sites, predicted_sites):
    """
    @param sites: path to sites.txt
    @param predicted_sites: path to predictedsites.txt
    @return norm: norm of the difference of sites and predicted sites
    @return count: number of overlapping sites
    """
    mySites = []
    myPredictedSites = []
    count = 0
    
    with open(sites, 'r') as sites_file:
        mySites = sites_file.readlines()
    with open(predicted_sites, 'r') as predicted_sites_file:
        myPredictedSites = predicted_sites_file.readlines()
    
    for i in range(len(mySites)):
        mySites[i] = int(mySites[i].strip())
        myPredictedSites[i] = int(myPredictedSites[i].strip())
        if mySites[i] == myPredictedSites[i]:
            count += 1
    
    sites_array = np.array(mySites)
    predicted_array = np.array(myPredictedSites)
    norm = la.norm(sites_array-predicted_array)
    
    return norm, count


directories = ['default', 'ICPC_1', 'ICPC_1.5', 'ML_6', 'ML_7', 'SC_5', 'SC_20'] 
dirname = os.path.dirname(__file__)
entropies = []
norm_list = []
count_list = []

for directory in directories:
    dir_entropy = []
    dir_norm = []
    dir_count = []
    for i in range(1, 11):
        predicted_motif_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/predictedmotif.txt')
        motif_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/motif.txt')
        sites_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/sites.txt')
        predicted_sites_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/predictedsites.txt')
        
        dir_entropy.append(relative_entropy(motif_file, predicted_motif_file))
        norm, count = sites_norm(sites_file, predicted_sites_file)
        dir_norm.append(norm)
        dir_count.append(count)
    entropies.append(sum(dir_entropy)/len(dir_entropy))
    norm_list.append(sum(dir_norm)/len(dir_norm))
    count_list.append(sum(dir_count)/len(dir_count))
        
print(entropies)
print(norm_list)
print(count_list)
       