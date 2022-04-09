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
    

sc = int(input("Enter sequence count: ")) #asks user for sequence count
sl = int(input("Enter sequence length: ")) #asks user for sequence length

mySequences = generate_sequences(sc, sl) #generate
print(mySequences) 


def generate_motif(ICPC, ML, SC):

    motif = []

    icpcDict = {1 : .8105, 1.5 : .9245, 2 : 1}
    randNum = np.random(0,1)
    p = icpcDict[ICPC]

    probA = (1-p)/3
    probC = (1-p)/3
    probG = (1-p)/3
    probT = (1-p)/3

    if (randNum < .25):
        probA = icpcDict[ICPC]
    elif (randNum < .55):
        probC = icpcDict[ICPC]
    elif (randNum < .75):
        probG = icpcDict[ICPC]
    else:
        probT = icpcDict[ICPC]

    for i in range(0, ML):
        seq = [0,0,0,0]
        for j in range(0, SC):
            randNum = np.random(0,1)
            if(randNum <= p):
                
            

    



    
    



