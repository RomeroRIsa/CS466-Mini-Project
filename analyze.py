from filecmp import dircmp
import os
import numpy as np
import numpy.linalg as la
import statistics


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

def gibbs_runtime(runtime):
    with open(runtime, 'r') as runtime_file:
        myRuntime = runtime_file.readlines()
    
    return float(myRuntime[0])


directories = ['default', 'ICPC_1', 'ICPC_1.5', 'ML_6', 'ML_7', 'SC_5', 'SC_20'] 
dirname = os.path.dirname(__file__)
entropies = []
norm_list = []
count_list = []
runtime_list = []
entropy_error = []
norm_error = []
runtime_error = []
count_error = []

for directory in directories:
    dir_entropy = []
    dir_norm = []
    dir_count = []
    dir_runtime = []
    dir_entropy_error = []
    dir_norm_error = []
    dir_runtime_error = []
    for i in range(1, 11):
        predicted_motif_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/predictedmotif.txt')
        motif_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/motif.txt')
        sites_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/sites.txt')
        predicted_sites_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/predictedsites.txt')
        runtime_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/runtime.txt')
        
        entropy = relative_entropy(motif_file, predicted_motif_file)
        dir_entropy.append(entropy)
        norm, count = sites_norm(sites_file, predicted_sites_file)
        dir_norm.append(norm)
        dir_count.append(count)
        runtime = gibbs_runtime(runtime_file)
        dir_runtime.append(runtime)
        
    entropies.append(sum(dir_entropy)/len(dir_entropy))
    norm_list.append(sum(dir_norm)/len(dir_norm))
    count_list.append(sum(dir_count)/len(dir_count))
    runtime_list.append(sum(dir_runtime)/len(dir_runtime))
    entropy_error.append(statistics.pstdev(dir_entropy)/np.sqrt(len(dir_entropy)))
    norm_error.append(statistics.pstdev(dir_norm)/np.sqrt(len(dir_norm)))
    runtime_error.append(statistics.pstdev(dir_runtime)/np.sqrt(len(dir_runtime)))
    count_error.append(statistics.pstdev(dir_count)/np.sqrt(len(dir_count)))
    
        
print("entropy:", entropies)
print("norm:", norm_list)
print("average overlapping sites:", count_list)
print("runtime:", runtime_list)
print("entropy error:", entropy_error)
print("norm error:", norm_error)
print("runtime error:", runtime_error)       
print("average overlapping sites error:", count_error)