# Run file for the Eulerian paths

import os
import time
import eulerian_path

min_k = 3
max_k = 31
sequence_count = 1000
run_times = {}

# Run all k values from min_k to max_k. Results are written to the contigs folder. Run times are recorded.
for k in range(min_k, max_k + 1):
    print(f'Running with k={k}: ', end='')
    start = time.time()
    eulerian_path.find_eulerian_path(sequence_count, k, f'../data/contigs/contig_k={k}.fna', f'k = {k}')
    end = time.time()
    run_times[k] = end - start
    print(end - start)

# Align each path created above. This makes runs a script that runs NeedlemanWunsch.java
for k in range(min_k, max_k + 1):
    print(f'Aligning k={k}...')
    os.system(f'java NeedlemanWunsch.java contig_k={k}.fna SARS-CoV-2_spike_protein.fna contig_alignment_k={k}.txt 1 -1 -2 0')
