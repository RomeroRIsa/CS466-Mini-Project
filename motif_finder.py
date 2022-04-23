from random import random
import numpy as np
import os
import math

def find_best_pos(sequence, motif, background, ML):
    best_pos = 0
    best_score = -np.inf
    nuc_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(len(sequence)-ML+1):
        score = 0
        for j in range(i, i+ML):
            nuc = nuc_dict[sequence[j]]
            if motif[j-i, nuc] == 0 or background[nuc] == 0:
                score = -np.inf
                break
            score += math.log(motif[j-i, nuc]/background[nuc])
        if score > best_score:
            best_score = score
            best_pos = i
    return best_pos


def gibbs_find_motif(sequence_strs, ML):
    sequence_list = np.array([[i for i in x] for x in sequence_strs])
    background = np.zeros((sequence_list.shape[0], 4))
    background[:,0] = np.count_nonzero(sequence_list=="A", axis=1)
    background[:,1] = np.count_nonzero(sequence_list=="C", axis=1)
    background[:,2] = np.count_nonzero(sequence_list=="G", axis=1)
    background[:,3] = np.count_nonzero(sequence_list=="T", axis=1)
    background = background/sequence_list.shape[1]
    sequence_pos = np.random.randint(0, sequence_list.shape[1]-ML+1, (sequence_list.shape[0]))
    prev_sequence_pos = np.full_like(sequence_pos, -1)
    motif = np.zeros((ML, 4))
    while not np.array_equal(prev_sequence_pos, sequence_pos):
        prev_sequence_pos = sequence_pos
        for i in range(sequence_list.shape[0]):
            idxs = ((np.arange(0, ML, dtype=np.int64) + np.zeros((sequence_list.shape[0],1), dtype=np.int64)) + sequence_pos[:, None])
            idxs = idxs + np.repeat(np.arange(0,sequence_list.shape[0]*sequence_list.shape[1], sequence_list.shape[1]), ML, axis=0).reshape(-1,ML)
            sequences = np.take(sequence_list, idxs)
            sequences = np.vstack((sequences[:i] ,sequences[i+1:])).T
            motif[:,0] = np.count_nonzero(sequences=="A", axis=1)
            motif[:,1] = np.count_nonzero(sequences=="C", axis=1)
            motif[:,2] = np.count_nonzero(sequences=="G", axis=1)
            motif[:,3] = np.count_nonzero(sequences=="T", axis=1)
            motif = motif/sequences.shape[1]
            sequence_pos[i] = find_best_pos(sequence_list[i], motif, background[i], ML)
            #print(motif)
    idxs = ((np.arange(0, ML, dtype=np.int64) + np.zeros((sequence_list.shape[0],1), dtype=np.int64)) + sequence_pos[:, None])
    idxs = idxs + np.repeat(np.arange(0,sequence_list.shape[0]*sequence_list.shape[1], sequence_list.shape[1]), ML, axis=0).reshape(-1,ML)
    sequences = np.take(sequence_list, idxs).T
    motif[:,0] = np.count_nonzero(sequences=="A", axis=1)
    motif[:,1] = np.count_nonzero(sequences=="C", axis=1)
    motif[:,2] = np.count_nonzero(sequences=="G", axis=1)
    motif[:,3] = np.count_nonzero(sequences=="T", axis=1)
    return motif/sequences.shape[1]
