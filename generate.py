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
            
    


