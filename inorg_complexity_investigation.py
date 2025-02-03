import signal
import time
import multiprocessing
import os
import tempfile
from cfg_util import *
from smiles_grammar import GCFG
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

# Function to process a single SMILES string
def process_single_smiles(smiles_data):
    smiles, idx = smiles_data
    time_limit = MAX_TIME

    # Set an alarm for the time limit
    signal.alarm(time_limit)

    temp_files = {}
    time_taken = None  # Initialize time_taken
    try:
        # Create temporary files for each time range dynamically
        for time_range in TIME_RANGES:
            temp_files[time_range] = tempfile.NamedTemporaryFile(delete=False, mode='w', prefix=f'batch_{idx}_range_{time_range[0]}_{time_range[1]}_')

        temp_files['failures'] = tempfile.NamedTemporaryFile(delete=False, mode='w', prefix=f'batch_{idx}_failures_')
        temp_files[f'>{MAX_TIME}'] = tempfile.NamedTemporaryFile(delete=False, mode='w', prefix=f'batch_{idx}_greater_than_{str(time_limit)}_')

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
            # Successful smiles taking less than 1 second are ignored (considered successful)
            return {'temp_files': [], 'time_taken': time_taken, 'successful': True}

        saved = False
        for time_range in TIME_RANGES:
            if time_range[0] <= time_taken < time_range[1]:
                temp_files[time_range].write(smiles + '\n')
                saved = True
                break

        # If time taken exceeds the maximum range of defined intervals
        if not saved and time_taken >= MAX_TIME:
            temp_files[f'>{MAX_TIME}'].write(smiles + '\n')
            # Close all temporary files
            for temp_file in temp_files.values():
                temp_file.close()
            return {
                'temp_files': [temp_files[f'>{MAX_TIME}'].name],
                'time_taken': time_taken,
                'successful': False
            }

        # Close all temporary files
        for temp_file in temp_files.values():
            temp_file.close()

        return {
            'temp_files': [temp_file.name for temp_file in temp_files.values() if temp_file.name],
            'time_taken': time_taken,
            'successful': True
        }

    except TimeoutException:
        # Handle the timeout exception and consider it a failure
        time_taken = time_limit  # Set time_taken to max limit for timeout
        temp_files[f'>{MAX_TIME}'].write(smiles + '\n')
        # Close all temporary files
        for temp_file in temp_files.values():
            temp_file.close()
        return {
            'temp_files': [temp_files[f'>{MAX_TIME}'].name],
            'time_taken': time_taken,
            'successful': False
        }
    except Exception:
        # Handle any other exceptions by writing to failures file
        temp_files['failures'].write(smiles + '\n')
        time_taken = None  # Set to None for unknown errors
        # Close all temporary files
        for temp_file in temp_files.values():
            temp_file.close()
        return {
            'temp_files': [temp_files['failures'].name],
            'time_taken': time_taken,
            'successful': False
        }
    finally:
        # Disable the alarm
        signal.alarm(0)

# Main function to execute multiprocessing
def main():
    start = time.time()
    smiles_file = 'smiles_inorganic.smi'
    total_processes = 8
    time_limit = MAX_TIME  # Set to the max time limit (30 seconds)

    # Read smiles from file and store them in a list
    with open(smiles_file, 'r') as f:
        smiles_list = f.readlines()
        smiles_list = [smiles.strip() for smiles in smiles_list]
        
    # For debugging purposes shorten the list
    smiles_list = smiles_list[:200000]

    # Create the data for multiprocessing (each SMILES with its index)
    smiles_data = [(smiles, idx) for idx, smiles in enumerate(smiles_list)]

    # Use a Pool for multiprocessing
    with multiprocessing.Pool(processes=total_processes) as pool:
        results = pool.map(process_single_smiles, smiles_data)

    # Aggregate results and merge temporary files into final files
    final_files = {time_range: f'final_range_{time_range[0]}_{time_range[1]}.txt' for time_range in TIME_RANGES}
    final_files['failures'] = 'final_failures.txt'
    final_files[f'>{MAX_TIME}'] = f'final_range_greater_than_{str(time_limit)}.txt'

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
                    if f'_greater_than_{str(time_limit)}_' in temp_file_name:
                        with open(temp_file_name, 'r') as temp_file:
                            final_file_handles[f'>{MAX_TIME}'].writelines(temp_file.readlines())
                    elif '_failures_' in temp_file_name or '_timeout_' in temp_file_name:
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
    max_batch_time = max(result['time_taken'] for result in results if result['time_taken'] is not None)
    total_processed = len(smiles_list)
    total_success = sum(1 for result in results if result['successful'])
    total_failures = smiles_count.get('failures', 0)
    total_greater_than_max = smiles_count.get(f'>{MAX_TIME}', 0)

    print("\nOverall Analysis:")
    print(f"{'Metric':<60}{'Count'}")
    print("-" * 70)
    print(f"{'Total Number of processed SMILES':<60}{total_processed}")
    print(f"{'Total Number of successful SMILES':<60}{total_success}")
    for key, count in smiles_count.items():
        if key not in ['failures', f'>{MAX_TIME}']:
            print(f"{f'Total Number of SMILES in {key}':<60}{count}")
    print(f"{'Total Number of SMILES in >MAX_TIME':<60}{total_greater_than_max}")
    print(f"{'Total Number of SMILES in failures':<60}{total_failures}")
    print(f"{'Max Time taken for any SMILES':<60}{max_batch_time:.2f} seconds")
    print(f"{'Total Time taken for program execution':<60}{time.time() - start:.2f} seconds")

if __name__ == '__main__':
    main()
