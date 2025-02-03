# Script to test the correct process of encoding a SMILES into a gene and decoding it back to a SMILES

from cfg_util import *
from smiles_grammar import GCFG
from GOs import mutation
import numpy as np
import time

# Start overall timing
overall_start = time.time()

# Initial SMILES input for testing
initial_smiles = 'O#C[Mo]123456(C=C71=C[Mo]189%1027(C#O)(C#O)c2c1c8c9c%102)(C#O)c1c3c4c5c61'
# Another example SMILES for variety
# initial_smiles = 'O=[Fe+]123N4=C(CN1(CC1=N2C(=CNC2CCCCC2)C=C1)CC1=N3C(=CNC2CCCCC2)C=C1)C=CC4=CNC1CCCCC1'
# Cl[Ru](C#O)([Si](Cl)(Cl)Cl)([P](C1CCCCC1)(C1CCCCC1)C1CCCCC1)[Si-](Cl)(Cl)Cl

smiles = initial_smiles

# --- Encoding ---
encoding_start = time.time()
encoded_smiles = encode(smiles)
encoding_end = time.time()
# print(f'Encoded smiles: {encoded_smiles}\n')


# --- From encoded SMILES to gene ---
gene_conversion_start = time.time()
gene = cfg_to_gene(encoded_smiles, max_len=-1)
gene_conversion_end = time.time()
# print(f'Gene: {gene}\n')

# --- From gene to decoded SMILES (CFG) ---
gene_to_cfg_start = time.time()
decoded_smiles = gene_to_cfg(gene)
gene_to_cfg_end = time.time()
# print(f'Decoded smiles: {decoded_smiles}\n')

# --- Decoding from CFG back to SMILES ---
cfg_to_smiles_start = time.time()
final_smiles = decode(decoded_smiles)
cfg_to_smiles_end = time.time()
# print(f'Final smiles: {final_smiles}\n')

# --- Total time taken ---
overall_end = time.time()

print(f'Time taken for encoding: {encoding_end - encoding_start:.2f} seconds\n')
print(f'Time taken for encoding to gene conversion: {gene_conversion_end - gene_conversion_start:.2f} seconds\n')
print(f'Time taken for gene to CFG conversion: {gene_to_cfg_end - gene_to_cfg_start:.2f} seconds\n')
print(f'Time taken for CFG to SMILES decoding: {cfg_to_smiles_end - cfg_to_smiles_start:.2f} seconds\n')

print(f'Total time taken: {overall_end - overall_start:.2f} seconds')
