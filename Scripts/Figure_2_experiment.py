"""
Script to generate the data to recreate Figure 2 in the manuscript: 
Benchmark spacer alignment time: comparing FuzzySearch vs. WFA (Wavefront Alignment)
with dictionary-based aligner and error-prone aligner.

This script generates synthetic spacer sequences with varying lengths and mismatch rates,
then measures the execution time of:
1. FuzzySearch (with different max_l_dist cutoffs)
2. WFA using a pre-constructed dictionary of aligners for quasi-local alignment (optimized)
3. WFA using a standard aligner

The results are saved to a CSV file.
"""

import csv
import os
import time
import random
import fuzzysearch
from pywfa import WavefrontAligner

firstRepeat = 'GAATTGAAAC'

def construct_aligners(pattern, min_text_length, max_text_length):
    aligner_dict = {}
    for length in range(min_text_length, max_text_length + 1, 10):
        aligner = WavefrontAligner(pattern, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, text_end_free = length, text_begin_free = length)
        aligner_dict[length] = aligner

    return aligner_dict

def dr1_sequence(sequence, MM_no, length, position):
    """
    Generates a synthetic read containing a repeat sequence with mismatches.
    
    Args:
        sequence (str): The repeat sequence to embed.
        MM_no (int): Number of mismatches to introduce into the repeat.
        length (int): Total length of the generated read.
        position (float): Relative position of the repeat in the read (0.0 to 1.0).
        
    Returns:
        str: The generated synthetic read.
    """
    start_position = int(((position * length)-(0.5 * len(sequence))))     
    num_nucleotides_before = start_position
    num_nucleotides_after = length - (num_nucleotides_before + len(sequence))  

    nucleotides_before = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_before))
    nucleotides_after = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_after))

    MM_positions = set()
    for _ in range(MM_no):
        while True:
            # We want to mutate the sequence itself
            pos = random.randint(0, len(sequence) - 1)
            if pos not in MM_positions:
                MM_positions.add(pos)
                break

        letter = random.choice('ATGC')
        sequence = sequence[:pos] + letter + sequence[pos + 1:]

    read = f"{nucleotides_before}{sequence}{nucleotides_after}"
    return read

def main():
    # --- 1. Construct Aligners ---
    # Pre-compute WFA aligners for different lengths to optimize performance
    # The range 10 to 150 covers the expected variability in read lengths (step size implicit in implementation)
    print("Constructing aligners...")
    dr1_aligner_dict = construct_aligners(firstRepeat, 10, 150)  
    
    # Create a standard aligner for comparison (representing less optimized or "error" case)
    dr1_aligner_error = WavefrontAligner(
        firstRepeat, 
        gap_opening=4, 
        gap_opening2=24, 
        gap_extension=4, 
        text_end_free=0, 
        text_begin_free=0
    )

    Method = []
    Length = []
    Time_ = []
    l_dist = []
    No_MM = []

    n_reads = 1000
    n_reps  = 5  # measure each length 5 times and average

    print("Starting benchmark loops...")

    for MM_no in [0, 1, 2, 3, 4, 5]: 
        print("MM: ", MM_no)

        # step size 5
        for length in range(50, 161, 5): 
            if length % 10 == 0:
                print("Length: ", length)

            # --- 2. Generate Sequences ---
            # Pre-generate the reads once per length (keeps comparability across methods/cutoffs)
            sequences = [dr1_sequence(firstRepeat, MM_no, length, 0.5) for _ in range(n_reads)]

            # --- 3. Benchmark FuzzySearch ---
            # score_cutoff 1..5
            for cutoff in [1, 2, 3, 4, 5]: 
                rep_times = []
                for _ in range(n_reps):
                    t0 = time.perf_counter()
                    for L in sequences:
                        _ = fuzzysearch.find_near_matches(firstRepeat, L, max_l_dist=cutoff)
                    t1 = time.perf_counter()
                    rep_times.append(t1 - t0)

                avg_time = sum(rep_times) / len(rep_times)

                Time_.append(avg_time)
                Method.append("fuzzy")
                Length.append(length)
                l_dist.append(cutoff)
                No_MM.append(MM_no)

            # --- 4. Benchmark WFA (Dictionary Aligner) ---
            # Quasi-local alignment
            rep_times = []
            for _ in range(n_reps):
                t0 = time.perf_counter()
                for L in sequences:
                    # Retrieve the correct aligner based on read length
                    aligner = dr1_aligner_dict[round(len(L), -1) - 10]
                    _ = aligner(L, clip_cigar=False)
                t1 = time.perf_counter()
                rep_times.append(t1 - t0)

            avg_time = sum(rep_times) / len(rep_times)

            Time_.append(avg_time)
            Method.append("wfa")
            Length.append(length)
            l_dist.append("-")
            No_MM.append(MM_no)

            # --- 5. Benchmark WFA (Standard Aligner) ---
            # Ends-free global alignment
            rep_times = []
            for _ in range(n_reps):
                t0 = time.perf_counter()
                for L in sequences:
                    _ = dr1_aligner_error(L, clip_cigar=False)
                t1 = time.perf_counter()
                rep_times.append(t1 - t0)

            avg_time = sum(rep_times) / len(rep_times)

            Time_.append(avg_time)
            Method.append("wfa")
            Length.append(length)
            l_dist.append("--")
            No_MM.append(MM_no)

    # --- 6. Write Results ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(os.path.dirname(script_dir), "Outputs", "Figure_2")
    os.makedirs(output_dir, exist_ok=True)

    csv_file = os.path.join(output_dir, "figure_2_experiment_data.csv")
    data = list(zip(Method, l_dist, Length, No_MM, Time_))

    with open(csv_file, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Method", "score_cutoff", "Length", "No_MM", "Time"])
        writer.writerows(data)

    print(f"CSV file '{csv_file}' created.")

if __name__ == "__main__":
    main()