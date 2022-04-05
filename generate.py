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


def generate_motif(ICPC, ML, seqList):
    
    motif = []

    for seq in seqList:
        for i in range(0, ML):
            aCount = 0
            cCount = 0
            gCount = 0
            tCount = 0
            if (seq[i] == "A"):
                aCount +=1
            if (seq[i] == "C"):
                cCount +=1
            if (seq[i] == "G"):
                gCount +=1
            if (seq[i] == "T"):
                tCount +=1

            total = aCount + cCount + gCount + tCount
            probA = aCount/total
            probC= cCount/total
            probG = gCount/total
            probT = tCount/total

            row = []
            row.append(probA)
            row.append(probC)
            row.append(probG)
            row.append(probT)
            motif.append(row)


