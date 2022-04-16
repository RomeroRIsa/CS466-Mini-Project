from random import random
import numpy as np
import os

def gibbs_find_motif(sequence, ML):
	sequence_list = np.array([[i for i in x] for x in sequence])
	background = np.zeros((sequence_list.shape[0], 4))
	background[:,0] = np.count_nonzero(sequence_list=="A", axis=1)
	background[:,1] = np.count_nonzero(sequence_list=="C", axis=1)
	background[:,2] = np.count_nonzero(sequence_list=="G", axis=1)
	background[:,3] = np.count_nonzero(sequence_list=="T", axis=1)
	background = background/sequence_list.shape[1]
	sequence_pos = np.random.randint(0, sequence_list.shape[1]-ML, (sequence_list.shape[0]))
	prev_sequence_pos = np.full_like(sequences, -1)
	while prev_sequence_pos != sequence_pos:
	    prev_sequence_pos = sequence_pos
	    for i in range(sequence_list.shape[0]):
	        idxs = (np.arange(0, ML, dtype=np.int64) + np.zeros((sequence_list.shape[0],1), dtype=np.int64)) + sequence_pos[:, None]
	        sequences = np.take(sequence_list, idxs)
	        sequences = np.vstack((sequences[:i] ,sequences[i+1:]))
	        motif = np.zeros((sequences.shape[0], 4))
	        motif[:,0] = np.count_nonzero(sequence_list=="A", axis=1)
	        motif[:,1] = np.count_nonzero(sequence_list=="C", axis=1)
	        motif[:,2] = np.count_nonzero(sequence_list=="G", axis=1)
	        motif[:,3] = np.count_nonzero(sequence_list=="T", axis=1)
	        motif = motif/sequences.shape[1]
	        
	        for 