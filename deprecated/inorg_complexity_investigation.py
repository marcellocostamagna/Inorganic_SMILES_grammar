import signal
import time
import multiprocessing
import os
import tempfile
from cfg_util import *
from smiles_grammar_test import GCFG
import numpy as np

# Define a timeout handler
class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException

# Set the signal handler for alarm
signal.signal(signal.SIGALRM, timeout_handler)

# Define time ranges for saving molecules (e.g., 1-5, 5-10, etc.)
TIME_RANGES = [(1, 5), (5, 10), (10, 15), (15, 20), (20, 25), (25, 30)]
MAX_TIME = 30  # Maximum time limit in seconds

# Function to process a batch of SMILES strings
def process_smiles_batch(batch_index, smiles_batch, time_limit, results_queue):
    start = time.time()

    # Create temporary files for each time range
    temp_files = {}
    for time_range in TIME_RANGES:
        temp_files[time_range] = tempfile.NamedTemporaryFile(delete=False, mode='w', prefix=f'batch_{batch_index}_range_{time_range[0]}_{time_range[1]}_')

    temp_files[f'>{MAX_TIME}'] = tempfile.NamedTemporaryFile(delete=False, mode='w', prefix=f'batch_{batch_index}_range_greater_than_{str(time_limit)}')
    temp_files['failures'] = tempfile.NamedTemporaryFile(delete=False, mode='w', prefix=f'batch_{batch_index}_failures_')

    successful_smiles = 0  # Count of successful SMILES in the batch

    for smiles in smiles_batch:
        # Set an alarm for the time limit
        signal.alarm(time_limit)
        
        try:
            # Start timing the processing of each SMILES
            start_partial = time.time()
            
            # Attempt to encode and decode the SMILES
            # Encoding
            encoded_smiles = encode(smiles)
            if encoded_smiles is None:
                raise ValueError("Encoding failed")

            # From encoded smiles to gene
            gene = cfg_to_gene(encoded_smiles, max_len=-1)

            # From gene to decoded smiles
            decoded_smiles = gene_to_cfg(gene)

            # Decoding (from decoded smiles back to smiles)
            final_smiles = decode(decoded_smiles)

            # Calculate the time taken
            time_taken = time.time() - start_partial

            # Categorize the smiles based on the time taken
            if time_taken < 1:
                # Ignore smiles that took less than 1 second and count them as successful
                successful_smiles += 1
                continue
            
            saved = False
            for time_range in TIME_RANGES:
                if time_range[0] <= time_taken < time_range[1]:
                    temp_files[time_range].write(smiles + '\n')
                    saved = True
                    break

            # If time taken exceeds the maximum range of defined intervals
            if not saved and time_taken >= MAX_TIME:
                temp_files[f'>{MAX_TIME}'].write(smiles + '\n')

            # Increment successful count if processed successfully
            if saved:
                successful_smiles += 1

        except TimeoutException:
            # Handle the timeout exception
            # Save SMILES that timed out in '>MAX_TIME' file
            temp_files[f'>{MAX_TIME}'].write(smiles + '\n')
        except Exception:
            # Handle any other exceptions by writing to failures file
            temp_files['failures'].write(smiles + '\n')
        finally:
            # Disable the alarm
            signal.alarm(0)

    # Close all temporary files
    for temp_file in temp_files.values():
        temp_file.close()

    # Record results in the queue for each batch
    results_queue.put({
        'batch_index': batch_index,
        'temp_files': [temp_file.name for temp_file in temp_files.values()],
        'time_taken': time.time() - start,
        'successful_smiles': successful_smiles,
    })

# Main function to execute multiprocessing
def main():
    smiles_file = 'smiles_inorganic.smi'
    total_processes = 8
    time_limit = MAX_TIME  # Set to the max time limit (30 seconds)

    # Read smiles from file and store them in a list
    with open(smiles_file, 'r') as f:
        smiles_list = f.readlines()
        smiles_list = [smiles.strip() for smiles in smiles_list]
        
    # for debugging purposes shorten the list
    smiles_list = smiles_list[:200]

    # Split the SMILES list into approximately equal chunks for each process
    chunk_size = len(smiles_list) // total_processes
    batches = [smiles_list[i * chunk_size: (i + 1) * chunk_size] for i in range(total_processes - 1)]
    batches.append(smiles_list[(total_processes - 1) * chunk_size:])  # The last batch takes the remainder

    processes = []
    results_queue = multiprocessing.Queue()

    # Start each process
    print(f'Starting {total_processes} processes...')
    for i, batch in enumerate(batches):
        p = multiprocessing.Process(target=process_smiles_batch, args=(i, batch, time_limit, results_queue))
        processes.append(p)
        p.start()

    # Collect results
    results = []
    for _ in processes:
        results.append(results_queue.get())

    # Wait for all processes to finish
    for p in processes:
        p.join()

    # Aggregate results and merge temporary files into final files
    final_files = {time_range: f'final_range_{time_range[0]}_{time_range[1]}.txt' for time_range in TIME_RANGES}
    final_files[f'>{MAX_TIME}'] = f'final_range_greater_than_{str(time_limit)}.txt'
    final_files['failures'] = 'final_failures.txt'

    # Open final files in append mode and merge temporary files into them
    final_file_handles = {key: open(file, 'a') for key, file in final_files.items()}
    
    try:
        # Go through all result temporary files and merge them
        for result in results:
            for temp_file_name in result['temp_files']:
                # Determine which range the temp file belongs to
                for time_range in TIME_RANGES:
                    if f'_range_{time_range[0]}_{time_range[1]}_' in temp_file_name:
                        with open(temp_file_name, 'r') as temp_file:
                            final_file_handles[time_range].writelines(temp_file.readlines())
                        break
                else:
                    if f'_range_greater_than_{str(time_limit)}' in temp_file_name:
                        with open(temp_file_name, 'r') as temp_file:
                            final_file_handles[f'>{MAX_TIME}'].writelines(temp_file.readlines())
                    elif '_failures_' in temp_file_name:
                        with open(temp_file_name, 'r') as temp_file:
                            final_file_handles['failures'].writelines(temp_file.readlines())

                # Delete the temporary file after merging
                os.remove(temp_file_name)
    finally:
        # Ensure all final files are properly closed
        for handle in final_file_handles.values():
            handle.close()

    # Count the number of SMILES in each final file
    smiles_count = {}
    for key, file_name in final_files.items():
        with open(file_name, 'r') as f: 
            smiles_count[key] = sum(1 for _ in f)

    # Print overall analysis in tabular format
    max_batch_time = max(result['time_taken'] for result in results)
    total_processed = len(smiles_list)
    total_success = sum(result['successful_smiles'] for result in results) 

    print("\nOverall Analysis:")
    print(f"{'Metric':<60}{'Count'}")
    print("-" * 70)
    print(f"{'Total Number of processed SMILES':<60}{total_processed}")
    print(f"{'Total Number of successful SMILES':<60}{total_success}")
    for key, count in smiles_count.items():
        print(f"{f'Total Number of SMILES in {key}':<60}{count}")
    print(f"{'Total Time taken for program execution':<60}{max_batch_time:.2f} seconds")

if __name__ == '__main__':
    main()
