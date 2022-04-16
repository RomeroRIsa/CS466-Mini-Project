from random import random
import numpy as np
import os


def generate_sequences(sequence_count, sequence_length):
    """
    @param sequence_count: number of sequences to be generated
    @param sequence_length: length of each seqeunce
    @return: sequences which is a list of strings of A, C, G, and Ts
    """
    sequences = [] #list of sequences
    sequence_array = np.random.rand(sequence_count, sequence_length) #np array of random number [0,1) of shape (sc, sl)
    
    for i in range(sequence_array.shape[0]):
        sequence_string = ""
        for j in range(sequence_array.shape[1]):
            if sequence_array[i][j] < 0.25:
                sequence_string += "A"
            elif sequence_array[i][j] < 0.5:
                sequence_string += "C"
                
            elif sequence_array[i][j] < 0.75:
                sequence_string += "G"
            else:
                sequence_string += "T"
        sequences.append(sequence_string)
            
    return sequences

def generate_binding_sites(sequence_count, sequence_length, motif_length):
    """
    @param sequence_count: number of sequences to be generated
    @param sequence_length: length of each seqeunce
    @param motif_length: length of each motif
    @return: np array of integers of size sequence count
    """
    
    return np.random.randint(0, sequence_length-motif_length-1, size=sequence_count)

def generate_motif(ICPC, ML, SC):
    """
    @param ICPC: information count per column
    @ML: motif length
    @SC: sequence count
    Used the algorithm given in the appendix of the mini project PDF
    @return: a list of lists with integers as the element
    """
    motif = []

    icpcDict = {1 : .8105, 1.5 : .9245, 2 : 1} #probabilities of each ICPC
    randNum = random()

    #select preferred nucleotide
    for i in range(0, ML):
        p = icpcDict[ICPC]
        otherP = (1-p)/3

        preffered = ""
        second = ""
        third = ""
        fourth = ""

        if (randNum < .25):
            preffered = "A"
            second = "C"
            third = "G"
            fourth = "T"

        elif (randNum < .50):
            preffered = "C"
            second = "A"
            third = "G"
            fourth = "T"
        elif (randNum < .75):
            preffered = "G"
            second = "A"
            third = "C"
            fourth = "T"
        else:
            preffered = "T"
            second = "A"
            third = "G"
            fourth = "C"

        #add 1 to nucleotide count based on probability
        seq = {"A" : 0, "G" : 0, "C" : 0, "T" : 0}
        for j in range(ML):
            randNum = random()
            if(randNum < p):
                seq[preffered] = p
            elif(randNum < p + otherP):
                seq[second] = otherP
            elif(randNum < p + 2*otherP):
                seq[third] = otherP
            else:
                seq[fourth] = otherP

        seqAdd = [seq["A"], seq["G"], seq["C"], seq["T"]]

        motif.append(seqAdd)

    return motif

def plant_sites(sequences, sites, motif):
    """
    @param sequences: a list of strings, each string is a sequence to be planted
    @sites: an np array of integers
    @motif: a list of lists with integers as elements (PWM)
    @return: a list of sequences with the motifs planted
    """
    nuc_dict = {0: 'A', 1:'C', 2:'G', 3:'T'}
    planted_sequences = []
    for i, (sequence, site) in enumerate(zip(sequences,sites)):
        seqlist = list(sequence)
        for i, nuc in enumerate(motif):
            seqlist[i+site] = nuc_dict[np.random.choice(np.arange(4), 1, p=nuc)[0]]
        planted_sequences.append(''.join(seqlist))
    return planted_sequences

def convert_motif(ML, motif, file):
    """
    @param ML: integer, motif length
    @param motif: a list of lists with integers as elements
    @param file: a string containing relative path to the file
    """
    with open(file, 'w') as f:
        f.write('>MOTIF1\t' + str(ML))
        for i in range(len(motif)):
            f.write('\n')
            for j in range(len(motif[i])):
                f.write(str(motif[i][j]) + '\t')
        f.write('\n<')

sl = 500 #default value for sequence length
sc_default = ['default', 'ICPC_1', 'ICPC_1.5', 'ML_6', 'ML_7', 'SC_5', 'SC_20'] #directories that will have SC = 10 (default)
#dict that maps each directory with tuple of parameters (ICPC, ML, SC)
directory_dict = {'default': (2, 8, 10),
                 'ICPC_1': (1, 8, 10),
                 'ICPC_1.5': (1.5, 8, 10),
                 'ML_6': (2, 6, 10),
                 'ML_7': (2, 7, 10),
                 'SC_5': (2, 8, 5),
                 'SC_20': (2, 8, 20) 
                 }
dirname = os.path.dirname(__file__) #path of generate.py

#generate 10 sequences and binding sites for directories in sc_default
for directory in sc_default:
    for i in range(1, 11):
        #join for relative path to file
        sequence_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/sequences.fasta') 
        site_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/sites.txt')
        motif_length_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/motiflength.txt')
        motif_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/motif.txt')
        
        #generate
        mySequences = generate_sequences(directory_dict[directory][2], sl)
        mySites = generate_binding_sites(directory_dict[directory][2], sl, directory_dict[directory][1])
        myMotifLength = directory_dict[directory][1]
        myMotif = generate_motif(directory_dict[directory][0], directory_dict[directory][1], directory_dict[directory][2])
        plantedSequences = plant_sites(mySequences, mySites, myMotif)
        
        #open files    
        fasta_file = open(sequence_file, "w")
        sites = open(site_file, "w") 
        motif_length = open(motif_length_file, "w")
        motif = open(motif_file, "w")
        
        #write into files
        for j in range(len(plantedSequences)):
            if j == len(plantedSequences)-1:
                fasta_file.write(">" + "Sequence " + str(j+1) + "\n" + plantedSequences[j]) 
                sites.write(str(mySites[j])) 
            else:
                fasta_file.write(">" + "Sequence " + str(j+1) + "\n" + plantedSequences[j] + "\n") 
                sites.write(str(mySites[j]) + "\n") 
        motif_length.write(str(myMotifLength))
        convert_motif(directory_dict[directory][1], myMotif, motif_file)
        
        #close files  
        fasta_file.close()
        sites.close()
        motif_length.close()








