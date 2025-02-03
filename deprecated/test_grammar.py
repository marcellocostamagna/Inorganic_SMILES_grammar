import signal
from cfg_util import *
from smiles_grammar_inorganic import GCFG
import numpy as np
import time

# Define a timeout handler
class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException

# Set the signal handler for alarm
signal.signal(signal.SIGALRM, timeout_handler)

start = time.time()

smiles_file = 'smiles.txt'

total_smiles = 1000

time_limit = 1 # Time limit for each SMILES in seconds

n_changed = 0
n_failed = 0
n_success = 0

failed = []
success = []

partial_times = []

# Read smiles from file and store them in a list
with open(smiles_file, 'r') as f:
    smiles_list = f.readlines()
    smiles_list = [smiles.strip() for smiles in smiles_list]
    
for smiles in smiles_list[:total_smiles]:
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
        
        # stop_partial = time.time()
        # partial_time = stop_partial - start_partial
        # partial_times.append(partial_time)
        
        # print(f'Partial time taken: {time.time() - start_partial} seconds')
    
    except TimeoutException:
        # Handle the timeout exception
        n_failed += 1
        print(f"Timeout occurred for SMILES: {smiles}")
    except Exception:
        # Handle any other exceptions silently
        n_failed += 1
        failed.append(smiles)
    finally:
        # Disable the alarm
        signal.alarm(0)
    
    
    
print(f'Number of successful conversions: {n_success}')
print(f'Number of failed conversions: {n_failed}')
print(f'Number of changed conversions: {n_changed}')

print(f'Efficiency: {n_success / total_smiles * 100}%')

# save failed and success SMILES
with open('failed.txt', 'w') as f:
    for smiles in failed:
        f.write(smiles + '\n')  
        
with open('success.txt', 'w') as f:
    for smiles in success:
        f.write(smiles + '\n')
        
# # print the mean of the partial times
# print(f'Mean partial time: {np.mean(partial_times)}')

print(f'Time taken: {time.time() - start} seconds')

