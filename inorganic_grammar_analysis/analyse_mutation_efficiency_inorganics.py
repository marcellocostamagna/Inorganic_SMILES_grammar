import signal
import time
import copy
import multiprocessing
import traceback
from ccdc.molecule import Molecule

# Switch to use the appropriate codebase
original_code = False
if original_code:
    from original_code.smiles_grammar import GCFG
else:    
    from smiles_grammar import GCFG

from cfg_util import *
from GOs import mutation

# Define a timeout handler
class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException

# Set the signal handler for alarm
signal.signal(signal.SIGALRM, timeout_handler)

# Function to process a batch of SMILES strings
def process_smiles_batch(smiles_batch, n_attempts, results_queue, encoding_time_limit=2):
    # Initializing counts for different types of success and failures
    n_success = 0
    n_unchanged = 0
    n_failed_empty = 0
    n_failed_real_error = 0
    n_valid = 0
    encoding_failures = 0
    encoding_timeout_failures = 0
    mutation_failures = 0
    decoding_failures = 0

    for smiles in smiles_batch:
        # Set an alarm for the time limit (encoding stage)
        signal.alarm(encoding_time_limit)
        try:
            # Encode the SMILES and convert to gene
            encoded_smiles = encode(smiles)
            gene = cfg_to_gene(encoded_smiles, max_len=-1)
        except TimeoutException:
            # Handle the timeout exception for encoding
            encoding_timeout_failures += 1
            # print(f"Encoding Timeout: {smiles}")
            continue
        except Exception as e:
            encoding_failures += 1
            # print(f"Encoding Failure: {smiles}")
            # print(traceback.format_exc())
            continue
        finally:
            # Disable the alarm
            signal.alarm(0)

        for _ in range(n_attempts):
            try:
                # MUTATION STEP
                mutated_gene = mutation(gene)
            except Exception as e:
                mutation_failures += 1
                print(f"Mutation Failure: Gene - {gene}")
                print(traceback.format_exc())
                continue

            try:
                # DECODING STEP
                mutated_decoded_smiles = gene_to_cfg(mutated_gene)
                new_smiles = decode(mutated_decoded_smiles)
            except Exception as e:
                decoding_failures += 1
                n_failed_real_error += 1
                print(f"Decoding Failure: Mutated Gene - {mutated_gene}")
                print(traceback.format_exc())
                continue

            # Tracking results
            if new_smiles == smiles:
                n_unchanged += 1
            elif new_smiles == '':
                n_failed_empty += 1
            elif new_smiles is None:
                decoding_failures += 1
                n_failed_real_error += 1
            elif new_smiles != smiles:
                n_success += 1

                # Validity check with CCDC 
                try:
                    mol = Molecule.from_string(new_smiles)
                    if mol is not None:
                        n_valid += 1
                except Exception as e:
                    n_failed_real_error += 1
                    # print(f"CCDC Validation Failure: SMILES - {new_smiles}")
                    # print(traceback.format_exc())
                    pass

    # Put the results in the queue
    results_queue.put({
        'n_success': n_success,
        'n_unchanged': n_unchanged,
        'n_failed_empty': n_failed_empty,
        'n_failed_real_error': n_failed_real_error,
        'n_valid': n_valid,
        'encoding_failures': encoding_failures,
        'encoding_timeout_failures': encoding_timeout_failures,
        'mutation_failures': mutation_failures,
        'decoding_failures': decoding_failures,
    })

# Main function to execute multiprocessing
def main():
    start_time = time.time()
    smiles_file = 'smiles_inorganic.smi'
    n_attempts = 100
    total_processes = 8  # Number of processes to use

    # Read smiles from file and store them in a list
    with open(smiles_file, 'r') as f:
        smiles_list = f.readlines()
        smiles_list = [smiles.strip() for smiles in smiles_list]
        
    # For debugging purposes shorten the list
    smiles_list = smiles_list[:100]

    # Filter out valid SMILES using CCDC
    valid_smiles = []
    for smiles in smiles_list:
        try:
            mol = Molecule.from_string(smiles)
            if mol is not None:
                valid_smiles.append(smiles)
        except:
            pass
            
    print(f'Number of starting valid smiles: {len(valid_smiles)}')

    # Create a multiprocessing queue to collect results
    results_queue = multiprocessing.Queue()

    # Start multiprocessing - assign SMILES evenly among processes
    smiles_per_process = len(valid_smiles) // total_processes
    processes = []

    for i in range(total_processes):
        start_idx = i * smiles_per_process
        end_idx = len(valid_smiles) if i == total_processes - 1 else (i + 1) * smiles_per_process
        batch = valid_smiles[start_idx:end_idx]
        p = multiprocessing.Process(target=process_smiles_batch, args=(batch, n_attempts, results_queue))
        processes.append(p)
        p.start()

    # Collect results
    results = []
    for _ in processes:
        results.append(results_queue.get())

    # Wait for all processes to complete
    for p in processes:
        p.join()

    # Aggregate results
    total_success = sum(result['n_success'] for result in results)
    total_unchanged = sum(result['n_unchanged'] for result in results)
    total_failed_empty = sum(result['n_failed_empty'] for result in results)
    total_failed_real_error = sum(result['n_failed_real_error'] for result in results)
    total_valid = sum(result['n_valid'] for result in results)
    total_encoding_failures = sum(result['encoding_failures'] for result in results)
    total_encoding_timeout_failures = sum(result['encoding_timeout_failures'] for result in results)
    total_mutation_failures = sum(result['mutation_failures'] for result in results)
    total_decoding_failures = sum(result['decoding_failures'] for result in results)

    total_processed = (len(valid_smiles) - total_encoding_timeout_failures) * n_attempts

    # Print overall analysis
    print("\nOverall Analysis:")
    print(f'Total number of attempts: {total_processed}')
    print(f'Number of successful mutations: {total_success}')
    print(f'Number of unchanged smiles: {total_unchanged}')
    print(f'Number of failed mutations due to empty SMILES: {total_failed_empty}')
    print(f'Number of failed mutations due to real errors: {total_failed_real_error}')
    print(f'Number of valid molecules: {total_valid}')
    print(f'Mutation efficiency (Syntactically correct): {total_success / total_processed * 100:.2f}%')
    print(f'Mutation efficiency (Chemically valid): {total_valid / total_processed * 100:.2f}%')
    print(f'Total Encoding Failures: {total_encoding_failures}')
    print(f'Total Encoding Timeout Failures: {total_encoding_timeout_failures}')
    print(f'Total Mutation Failures: {total_mutation_failures}')
    print(f'Total Decoding Failures: {total_decoding_failures}')
    print(f'Time taken: {time.time() - start_time:.2f} seconds')

if __name__ == '__main__':
    main()