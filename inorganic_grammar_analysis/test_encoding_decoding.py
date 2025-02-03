# Script to test the correct process of encoding a smiles into a gene and decoding it back to a smiles

from cfg_util import *
from smiles_grammar_inorganic import GCFG
from GOs import mutation
import numpy as np
import time

start = time.time()
#'CC(=O)OCC[N+](C)(C)C'
# initial_smiles = 'O=[Fe+]123N4=C(CN1(CC1=N2C(=CNC2CCCCC2)C=C1)CC1=N3C(=CNC2CCCCC2)C=C1)C=CC4=CNC1CCCCC1' #Does not finish
initial_smiles = 'O#C[Mo]123456(C=C71=C[Mo]189%1027(C#O)(C#O)c2c1c8c9c%102)(C#O)c1c3c4c5c61' #4.5 seconds
# initial_smiles = 'O#C[Fe]1234(C#O)(C#CC5=CC=CC6=NSN=C56)c5c1c2c3c45' #3.5
# initial_smiles = 'O=[Cl](=O)(=O)O[Zn+]123N4=C(CN1(CC1=N2C(=CNC2CCCCC2)C=C1)CC1=N3C(=CNC2CCCCC2)C=C1)C=CC4=CNC1CCCCC1' #DOesnot finish
# initial_smiles = 'Cl[Ru](C#O)([Si](Cl)(Cl)Cl)([P](C1CCCCC1)(C1CCCCC1)C1CCCCC1)[Si-](Cl)(Cl)Cl' # 0.05 seconds
# initial_smiles = 'c1cc(ccc1C=Cc12c3c4c5c1[Fe]16782345c2c1c6c7c82)C=Cc12c3c4c5c1[Fe]16782345c2c1c6c7c82' # 15 seconds
# initial_smiles = 'c1cc(ccc1C=Cc12c3c4c5c1[Fe]16782345c2c1c6c7c82)' # 0.12 seconds
# initial_smiles = 'c1cc(ccc1C=Cc12c3c4c5c1[Fe]16782345c2c1c6c7c82)C=Cc12c3c4c5c1[Fe]2345' # 3.6 seconds
# initial_smiles = 'c1cc(ccc1C=Cc12c3c4c5c1[Fe]16782345c2c1c6c7c82)C=Cc12c3c4c5c1[Fe]2345CCCCC' # 14 seconds
# initial_smiles = 'c1cc(ccc1C=Cc12c3c4c5c1[Fe]16782345c2c1c6c7c82)C=Cc12c3c4c5c1[Fe]2345CCCCC(C#O)'
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
# print(f'Original smiles: {smiles}')
print(f'Final smiles:    {final_smiles}\n')

print(f'Time taken: {time.time() - start} seconds')













