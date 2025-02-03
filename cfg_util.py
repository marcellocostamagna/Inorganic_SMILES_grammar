import nltk

import nltk.parse.chart
import numpy as np

from smiles_grammar import GCFG


def get_smiles_tokenizer(cfg):
    long_tokens = [a for a in cfg._lexical_index.keys() if len(a) > 1]
    # there are currently 6 double letter entities in the grammar (7 with a new metal)
    # these are their replacement, with no particular meaning
    # they need to be ascii and not part of the SMILES symbol vocabulary
    # replacements = ['!', '?', '.', ',', ';', '$', '_'] #(one symbol added)
    
    replacements = [
    '!', '?', '$', '&', '*', '~', '_', ';', '.', ',',
    '`', '|', '<', '>', '{', '}', '§', '^', 'a', 'A',
    'z', 'Z', 'm', 'M', 'd', 'D', 'e', 'E', 'g', 'G',
    'j', 'J', 'l', 'L', 'q', 'Q', 'r', 'R', 't', 'T',
    'x', 'X', '"', '´', '˚', 'å', 'Å', 'ø', 'Ø', '¨',
    'ä', 'Ä', 'ö', 'Ö', 'ü', 'Ü', 'á', 'Á', 'é', 'É',
    'í', 'Í', 'ó', 'Ó', 'ú', 'Ú', 'à', 'À', 'è', 'È',
    'ì', 'Ì', 'ò', 'Ò', 'ù', 'Ù', 'â', 'Â', 'ê', 
    ]
    
    
    # print(f'Long tokens: {len(long_tokens)}')
    # # print(f'Long tokens: {long_tokens}')
    # print(f'Replacements: {len(replacements)}')
    assert len(long_tokens) == len(replacements)
    for token in replacements:
        assert token not in cfg._lexical_index

    def tokenize(smiles):
        for i, token in enumerate(long_tokens):
            smiles = smiles.replace(token, replacements[i])
        tokens = []
        for token in smiles:
            try:
                ix = replacements.index(token)
                tokens.append(long_tokens[ix])
            except Exception:
                tokens.append(token)
        return tokens
    
    return tokenize


def encode(smiles):
    tokenize = get_smiles_tokenizer(GCFG)
    tokens = tokenize(smiles)
    parser = nltk.ChartParser(GCFG)
    try:
        parse_tree = parser.parse(tokens).__next__()
    except StopIteration:
        # print(f'Failed to parse {smiles}')
        return None
    productions_seq = parse_tree.productions()
    # print(productions_seq)
    # print(parse_tree.productions())
    # print(f'Length of productions_seq: {len(productions_seq)}')
    productions = GCFG.productions()
    prod_map = {}
    for ix, prod in enumerate(productions):
        prod_map[prod] = ix
    indices = np.array([prod_map[prod] for prod in productions_seq], dtype=int)  
    return indices


def prods_to_eq(prods):
    seq = [prods[0].lhs()]
    for prod in prods:
        if str(prod.lhs()) == 'Nothing':
            break
        for ix, s in enumerate(seq):
            if s == prod.lhs():
                seq = seq[:ix] + list(prod.rhs()) + seq[ix + 1:]
                break
    try:
        return ''.join(seq)
    except Exception:
        return ''

def decode(rule):
    productions = GCFG.productions()
    prod_seq = [productions[i] for i in rule]
    # print(prod_seq)
    # print(f'Length of prod_seq: {len(prod_seq)}')
    return prods_to_eq(prod_seq)

def cfg_to_gene(prod_rules, max_len=-1):
    gene = []
    for r in prod_rules:
        lhs = GCFG.productions()[r].lhs()
        possible_rules = [idx for idx, rule in enumerate(GCFG.productions())
                          if rule.lhs() == lhs]
        gene.append(possible_rules.index(r))
    if max_len > 0:
        if len(gene) > max_len:
            gene = gene[:max_len]
        else:
            gene = gene + [np.random.randint(0, 256)
                           for _ in range(max_len - len(gene))]
    return gene


def gene_to_cfg(gene):
    prod_rules = []
    stack = [GCFG.productions()[0].lhs()]
    for g in gene:
        try:
            lhs = stack.pop()
        except Exception:
            break
        possible_rules = [idx for idx, rule in enumerate(GCFG.productions())
                          if rule.lhs() == lhs]
        rule = possible_rules[g % len(possible_rules)]
        prod_rules.append(rule)
        rhs = filter(lambda a: (type(a) == nltk.grammar.Nonterminal) and (str(a) != 'None'),
                     GCFG.productions()[rule].rhs())
        # rhs_list = list(rhs)
        # rhs_list_reversed = rhs_list[::-1]
        # stack.extend(rhs_list_reversed)
        stack.extend(list(rhs)[::-1])   
    return prod_rules
