import time
from cfg_util import *
from smiles_grammar_inorganic import GCFG
import numpy as np

start_general = time.time()

def process_smiles(smiles):
    start = time.time()
    # Encoding
    encoded_smiles = encode(smiles)

    # From encoded smiles to gene
    gene = cfg_to_gene(encoded_smiles, max_len=-1)

    # From gene to decoded smiles
    decoded_smiles = gene_to_cfg(gene)

    # Decoding (from decoded smiles back to smiles)
    final_smiles = decode(decoded_smiles)
    
    time_taken = time.time() - start
    
    return time_taken

range_fle = 'final_range_5_10.txt'
output_less = 'less_time.txt'
output_more = 'more_time.txt'

range = (5, 10)

# read smiles from file and make a list
with open(range_fle, 'r') as f:
    smiles = f.readlines()
    smiles_list = [smile.strip() for smile in smiles]
    
print(f'Number of smiles: {len(smiles_list)}')
correct = 0 

smiles_list = smiles_list[:10]
print(f'Number of smiles: {len(smiles_list)}')
    
for smile in smiles_list:
    time_taken = process_smiles(smile)
    # check if time is less, more or equal to the range
    if time_taken < range[0]:
        # save the smile to the less time file
        with open(output_less, 'a') as f:
            f.write(f'{smile} - {time_taken}\n')
    elif time_taken > range[1]:
        # save the smile to the more time file
        with open(output_more, 'a') as f:
            f.write(f'{smiles} - {time_taken}\n')
    else:
        correct+=1
    
print(f'Number of correct time estimates: {correct}')
print(f'Total time taken: {time.time() - start_general} seconds')
