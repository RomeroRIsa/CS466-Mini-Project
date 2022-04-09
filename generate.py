from random import random
import numpy as np
import os


def generate_sequences(sequence_count, sequence_length):
    """
    @param sequence_count: number of sequences to be generated
    @param sequence_length: length of each seqeunce
    """
    sequences = [] #list of sequences
    sequence_array = np.random.rand(sequence_count, sequence_length) #np array of random number [0,1)
    
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
    """
    
    return np.random.randint(0, sequence_length-motif_length-1, size=sequence_count)
    

sl = 500 #default value for sequence length
sc = [5, 10, 20] #possible sequence counts
ml = [6, 7, 8] #possible motif lengths
sc_default = ['default', 'ICPC_1', 'ICPC_1.5', 'ML_6', 'ML_7'] #directories that will have SC = 10 (default)
dirname = os.path.dirname(__file__) #path of generate.py

#generate 10 sequences and binding sites for directories in sc_default
for directory in sc_default:
    for i in range(1, 11):
        #join for relative path to file
        sequence_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/sequences.fasta') 
        site_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/sites.txt')
        motif_length_file = os.path.join(dirname, 'dataset/'+ directory + '/' + str(i) + '/motiflength.txt')
        
        mySequences = generate_sequences(sc[1], sl) #generate sequence
        
        if directory == 'ML_6':
            mySites = generate_binding_sites(sc[1], sl, ml[0])
            myMotifLength = ml[0]
        elif directory == 'ML_7':
            mySites = generate_binding_sites(sc[1], sl, ml[1])
            myMotifLength= ml[1]
        else:
            mySites = generate_binding_sites(sc[1], sl, ml[2])
            myMotifLength = ml[2]
        
        #open files    
        fasta_file = open(sequence_file, "w")
        sites = open(site_file, "w") 
        motif_length = open(motif_length_file, "w")
        
        for j in range(len(mySequences)):
            fasta_file.write(">" + "Sequence " + str(j+1) + "\n" + mySequences[j] + "\n") #write to fasta file
            sites.write(str(mySites[j]) + "\n") #write to sites.txt
        motif_length.write(str(myMotifLength))
          
        #close files  
        fasta_file.close()
        sites.close()

#generate 10 sequences and binding sites for SC = 5
for i in range(1, 11):
    #join for relative path to file
    sequence_file = os.path.join(dirname, 'dataset/SC_5/' + str(i) + '/sequences.fasta')
    site_file = os.path.join(dirname, 'dataset/SC_5/' + str(i) + '/sites.txt')
    motif_length_file = motif_length_file = os.path.join(dirname, 'dataset/SC_5/' + str(i) + '/motiflength.txt')
    
    mySequences = generate_sequences(sc[0], sl) #generate sequence
    mySites = generate_binding_sites(sc[0], sl, ml[2])
    myMotifLength = ml[2]
    
    
    #open files
    fasta_file = open(sequence_file, "w")
    sites = open(site_file, "w") 
    motif_length = open(motif_length_file, "w") 
    
    for j in range(len(mySequences)):
        fasta_file.write(">" + "Sequence " + str(j+1) + "\n" + mySequences[j] + "\n") #write to fasta file
        sites.write(str(mySites[j]) + "\n") #write to sites.txt
    motif_length.write(str(myMotifLength))
    
    #close files
    fasta_file.close()
    sites.close()
    motif_length.close()

#generate 10 sequences and binding sites for SC = 20
for i in range(1, 11):
    #join for relative path to file
    sequence_file = os.path.join(dirname, 'dataset/SC_20/' + str(i) + '/sequences.fasta')
    site_file = os.path.join(dirname, 'dataset/SC_20/' + str(i) + '/sites.txt')
    motif_length_file = os.path.join(dirname, 'dataset/SC_20/' + str(i) + '/motiflength.txt')
    
    mySequences = generate_sequences(sc[2], sl) #generate sequence
    mySites = generate_binding_sites(sc[2], sl, ml[2])
    myMotifLength = ml[2]
    
    #open files
    fasta_file = open(sequence_file, "w") 
    sites = open(site_file, "w")
    motif_length = open(motif_length_file, "w") 
    
    for j in range(len(mySequences)):
        fasta_file.write(">" + "Sequence " + str(j+1) + "\n" + mySequences[j] + "\n") #write to fasta file
        sites.write(str(mySites[j]) + "\n")
    motif_length.write(str(myMotifLength))
      
    #close files  
    fasta_file.close()
    sites.close()
    motif_length.close()


def plant_sites(sequences, sites, motif):
    ML = np.sum(motif[0])
    motif = motif/ML
    nuc_dict = {0: 'A', 1:'C', 2:'G', 3:'T'}
    planted_sequences = []
    for i, (sequence, site) in enumerate(zip(sequences,sites)):
        seqlist = list(sequence)
        for i, nuc in enumerate(motif):
            seqlist[i+site] = nuc_dict[np.random.choice(np.arange(4), 1, p=nuc)[0]]
        planted_sequences.append(''.join(seqlist))
    return planted_sequences


#sc = int(input("Enter sequence count: ")) #asks user for sequence count
#sl = int(input("Enter sequence length: ")) #asks user for sequence length

#mySequences = generate_sequences(sc, sl) #generate
#print(mySequences) 

ML = 5
SC = 20
ICPC = 1

def generate_motif(ICPC, ML, SC):
    
    motif = []

    icpcDict = {1 : .8105, 1.5 : .9245, 2 : 1}
    randNum = random()

    for i in range(0, ML):
        p = icpcDict[ICPC]
        otherP = (1-p)/3

        probA = (1-p)/3
        probC = (1-p)/3
        probG = (1-p)/3
        probT = (1-p)/3

        preffered = ""
        second = ""
        third = ""
        fourth = ""

        if (randNum < .25):
            probA = icpcDict[ICPC]
            preffered = "A"
            second = "C"
            third = "G"
            fourth = "T"

        elif (randNum < .50):
            probC = icpcDict[ICPC]
            preffered = "C"
            second = "A"
            third = "G"
            fourth = "T"
        elif (randNum < .75):
            probG = icpcDict[ICPC]
            preffered = "G"
            second = "A"
            third = "C"
            fourth = "T"
        else:
            probT = icpcDict[ICPC]
            preffered = "T"
            second = "A"
            third = "G"
            fourth = "C"


        seq = {"A" : 0, "G" : 0, "C" : 0, "T" : 0}
        for j in range(0, SC):
            randNum = random()
            if(randNum < p):
                seq[preffered] += 1
            elif(randNum < p + otherP):
                seq[second] += 1
            elif(randNum < p + 2*otherP):
                seq[third] += 1
            else:
                seq[fourth] += 1


        seqAdd = [seq["A"], seq["G"], seq["C"], seq["T"]]

        motif.append(seqAdd)

    return motif

motif = (generate_motif(ICPC, ML, SC))

def convert_motif(ML, motif):
    with open('motif.txt', 'w') as f:
        f.write('<MOTIF1\t' + str(ML))
        for i in range(len(motif)):
            f.write('\n')
            for j in range(len(motif[i])):
                f.write(str(motif[i][j]) + '\t')
        f.write('\n<')

convert_motif(ML, motif)


print(generate_motif(1, 6, 20))





