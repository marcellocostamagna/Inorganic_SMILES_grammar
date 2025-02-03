# Description: This file contains the functions for the genetic operators used in the genetic algorithm.
from rdkit import rdBase
rdBase.DisableLog('rdApp.*')
import numpy as np
import nltk
import copy
from smiles_grammar_inorganic import GCFG
from cfg_util import *

def mutation(gene):
    idx = np.random.choice(len(gene))
    gene_mutant = copy.deepcopy(gene)
    gene_mutant[idx] = np.random.randint(0, 256)
    return gene_mutant


def deduplicate(population):
    unique_smiles = set()
    unique_population = []
    for item in population:
        score, smiles, gene = item
        if smiles not in unique_smiles:
            unique_population.append(item)
        unique_smiles.add(smiles)
    return unique_population


def mutate(p_gene):
    c_gene = mutation(p_gene)
    try:
        # Decode the mutated gene into SMILES directly without RDKit canonicalization
        # (Devation from original Guacamol code)
        c_smiles = decode(gene_to_cfg(c_gene))
        
    except Exception as e:
        # Handle any decoding errors gracefully
        print(f'The mutation resulted in an invalid molecule: {e}')
        c_smiles = ''
        
    return c_smiles
