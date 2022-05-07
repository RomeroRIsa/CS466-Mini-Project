# CS466-Mini-Project
Intro to Bioinformatics Mini Project

Link to Google Docs: [Paper](https://docs.google.com/document/d/1-McSf4ZT1TPNfn0_Er0sGCpLtRS3tB_HOOIGxq2OvKc/edit#)

## generate.py
Generates 7 sets of 10 datasets, where each set has different parameters. 
- default: ICPC = 2, ML = 8, SC = 10, SL = 500
- ICPC_1: ICPC = 1, ML = 8, SC = 10, SL = 500
- ICPC_1.5: ICPC = 1.5, ML = 8, SC = 10, SL = 500
- ML_6: ICPC = 2, ML = 6, SC = 10, SL = 500
- ML_7: ICPC = 2, ML = 7, SC = 10, SL = 500
- SC_5: ICPC = 2, ML = 8, SC = 5, SL = 500
- SC_20: ICPC = 2, ML = 8, SC = 20, SL = 500

To run generate.py, run the following command in the terminal:
```python generate.py```

## motif_finder.py
Use the Gibbs sampling algorithm on the benchmark dataset to get "predictedmotif.txt" and "predictedsites.txt" for every dataset

To run motif_finder.py on the benchmark in the "dataset" directory, run the following command in the terminal:
```python motif_finder.py```

## analyze.py
Analyze the performance of the motif finding algorithm implemented in motif_finder.py by comparing "predictedmotif.txt" and "predictedsites.txt" with data generated in generate.py. Outputs the numbers for our measures of performance. 

To run analyze.py, run the following command in the terminal:
```python analyze.py```
