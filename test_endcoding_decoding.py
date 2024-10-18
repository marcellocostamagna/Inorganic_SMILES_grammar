# Script to test the correct process of encoding a smiles into a gene and decoding it back to a smiles

from cfg_util import *
from smiles_grammar_test import GCFG
from GOs import mutation
import numpy as np
import time

start = time.time()
#'CC(=O)OCC[N+](C)(C)C'
initial_smiles = 'CC(c1c(CC)cc(C=O)cc1)(CC(CO)CC)'
#Cl[Ru](C#O)([Si](Cl)(Cl)Cl)([P](C1CCCCC1)(C1CCCCC1)C1CCCCC1)[Si-](Cl)(Cl)Cl
smiles = initial_smiles

# Encoding
encoded_smiles = encode(smiles)
print(f'Encoded smiles: {encoded_smiles}\n')

# From encoded smiles to gene
gene = cfg_to_gene(encoded_smiles, max_len=-1)
print(f'Gene: {gene}\n')

# From gene to decoded smiles
decoded_smiles = gene_to_cfg(gene)
print(f'Decoded smiles: {decoded_smiles}\n')

# Decoding (from decoded smiles back to smiles)
final_smiles = decode(decoded_smiles)
print(f'Original smiles: {smiles}')
print(f'Final smiles:    {final_smiles}\n')

print(f'Time taken: {time.time() - start} seconds')













