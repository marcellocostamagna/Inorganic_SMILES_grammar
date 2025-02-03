import time
import os
import multiprocessing
from cfg_util import *
from smiles_grammar_inorganic import GCFG
import numpy as np

start_general = time.time()

# Define the range file and output files
range_file = 'final_range_1_5.txt'
output_less = 'less_time_parallel.txt'
output_more = 'more_time_parallel.txt'
range_ = (1, 5)

def process_smiles(smiles):
    start = time.time()
    try:
        # Encoding
        encoded_smiles = encode(smiles)

        # From encoded smiles to gene
        gene = cfg_to_gene(encoded_smiles, max_len=-1)

        # From gene to decoded smiles
        decoded_smiles = gene_to_cfg(gene)

        # Decoding (from decoded smiles back to smiles)
        final_smiles = decode(decoded_smiles)

    except Exception as e:
        print(f"Error processing SMILES: {smiles}, Error: {e}")
        return smiles, None  # Return None if processing fails

    time_taken = time.time() - start
    return smiles, time_taken

def process_smiles_worker(smiles_list, time_range, result_queue):
    results = []
    for smiles in smiles_list:
        processed_smiles, time_taken = process_smiles(smiles)
        if time_taken is not None:
            if time_taken < time_range[0]:
                results.append((processed_smiles, 'less'))
            elif time_taken > time_range[1]:
                results.append((processed_smiles, 'more'))
            else:
                results.append((processed_smiles, 'correct'))
    result_queue.put(results)

def main():
    # Read smiles from file and create a list
    with open(range_file, 'r') as f:
        smiles = f.readlines()
        smiles_list = [smile.strip() for smile in smiles]

    print(f'Number of SMILES: {len(smiles_list)}')
    
    smiles_list = smiles_list[:10]
    print(f'Number of SMILES: {len(smiles_list)}')

    # Use multiprocessing
    total_processes = multiprocessing.cpu_count()  # Utilize all available CPU cores
    chunk_size = len(smiles_list) // total_processes

    # Split the SMILES list into chunks for parallel processing
    chunks = [smiles_list[i:i + chunk_size] for i in range(0, len(smiles_list), chunk_size)]

    result_queue = multiprocessing.Queue()
    processes = []

    # Start worker processes
    for chunk in chunks:
        p = multiprocessing.Process(target=process_smiles_worker, args=(chunk, range_, result_queue))
        processes.append(p)
        p.start()

    # Collect results
    all_results = []
    for _ in processes:
        all_results.extend(result_queue.get())

    # Wait for all worker processes to finish
    for p in processes:
        p.join()

    # Write mismatched SMILES to respective files
    correct_count = 0
    with open(output_less, 'a') as less_file, open(output_more, 'a') as more_file:
        for smiles, category in all_results:
            if category == 'less':
                less_file.write(f"{smiles}\n")
            elif category == 'more':
                more_file.write(f"{smiles}\n")
            elif category == 'correct':
                correct_count += 1

    print(f'Number of correct time estimates: {correct_count}')
    print(f'Total time taken: {time.time() - start_general} seconds')

if __name__ == '__main__':
    main()
