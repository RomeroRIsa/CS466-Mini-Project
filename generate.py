from random import random
import numpy as np

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


    




            

    



    
    



