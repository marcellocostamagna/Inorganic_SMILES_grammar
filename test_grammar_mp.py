import signal
import time
import multiprocessing
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

# Function to process a batch of SMILES strings
def process_smiles_batch(batch_index, smiles_batch, time_limit, results_queue):
    start = time.time()
    n_changed, n_failed, n_success, n_timeout = 0, 0, 0, 0
    failed, success, timeout_smiles = [], [], []

    for smiles in smiles_batch:
        start_partial = time.time()
        
        # Set an alarm for the time limit
        signal.alarm(time_limit)
        
        try:
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

            # Check if the final smiles is the same as the original smiles
            if smiles == final_smiles:
                n_success += 1
                success.append(smiles)
            else:
                n_changed += 1

        except TimeoutException:
            # Handle the timeout exception
            n_timeout += 1
            timeout_smiles.append(smiles)
            # print(f"Timeout occurred for SMILES: {smiles}")
        except Exception:
            # Handle any other exceptions silently
            n_failed += 1
            failed.append(smiles)
        finally:
            # Disable the alarm
            signal.alarm(0)

    # Record results in the queue for each batch
    results_queue.put({
        'batch_index': batch_index,
        'n_success': n_success,
        'n_failed': n_failed,
        'n_changed': n_changed,
        'n_timeout': n_timeout,
        'failed': failed,
        'success': success,
        'timeout': timeout_smiles,
        'time_taken': time.time() - start
    })

# Main function to execute multiprocessing
def main():
    smiles_file = 'smiles.txt'
    total_processes = 7
    batch_size = 1000
    time_limit = 1  # Time limit for each SMILES in seconds
    start_index = 200000  # Starting index for processing the SMILES list

    # Read smiles from file and store them in a list
    with open(smiles_file, 'r') as f:
        smiles_list = f.readlines()
        smiles_list = [smiles.strip() for smiles in smiles_list]

    # Determine the number of batches and adjust if necessary
    remaining_smiles = smiles_list[start_index:]
    total_batches = min(total_processes, len(remaining_smiles) // batch_size)

    batches = [remaining_smiles[i * batch_size: (i + 1) * batch_size] for i in range(total_batches)]
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
        result = results_queue.get()
        results.append(result)

    # Wait for all processes to finish
    for p in processes:
        p.join()

    # Aggregate results and save them
    total_success = sum(result['n_success'] for result in results)
    total_failed = sum(result['n_failed'] for result in results)
    total_changed = sum(result['n_changed'] for result in results)
    total_timeout = sum(result['n_timeout'] for result in results)
    total_processed = total_success + total_failed + total_changed + total_timeout
    max_batch_time = max(result['time_taken'] for result in results)

    for result in results:
        batch_index = result['batch_index']
        
        # Save failed SMILES for each batch
        with open(f'failed_batch_{batch_index}.txt', 'w') as f:
            for smiles in result['failed']:
                f.write(smiles + '\n')
        
        # Save successful SMILES for each batch
        with open(f'success_batch_{batch_index}.txt', 'w') as f:
            for smiles in result['success']:
                f.write(smiles + '\n')
        
        # Save timed out SMILES for each batch
        with open(f'timeout_batch_{batch_index}.txt', 'w') as f:
            for smiles in result['timeout']:
                f.write(smiles + '\n')

        # Print batch summary
        print(f'\nBatch {batch_index} Summary:')
        print(f'  Successful conversions: {result["n_success"]}')
        print(f'  Failed conversions: {result["n_failed"]}')
        print(f'  Changed conversions: {result["n_changed"]}')
        print(f'  Timed out conversions: {result["n_timeout"]}')
        print(f'  Time taken for batch {batch_index}: {result["time_taken"]} seconds')
        
    # Print overall analysis
    print("\nOverall Analysis:")
    print(f'Total Number of successful conversions: {total_success}')
    print(f'Total Number of failed conversions: {total_failed}')
    print(f'Total Number of changed conversions: {total_changed}')
    print(f'Total Number of timed out conversions: {total_timeout}')
    print(f'Total Number of processed molecules: {total_processed}')
    print(f'Efficiency: {(total_success / total_processed) * 100}%')
    print(f'Total Time taken for program execution (longest batch time): {max_batch_time} seconds')
    

if __name__ == '__main__':
    main()
