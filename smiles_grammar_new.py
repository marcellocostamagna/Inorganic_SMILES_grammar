import nltk

# smiles grammar
gram = """
smiles -> chain
chain -> branched_atom
chain -> chain branched_atom
chain -> chain bond branched_atom

branched_atom -> atom
branched_atom -> atom RB
branched_atom -> atom BB
branched_atom -> atom RB BB

# Defining atom types
atom -> bracket_atom
atom -> aliphatic_organic
atom -> aromatic_organic

# Treating transition metals as generic bracket atoms
# Making [metal_symbol] a single recognized entity.
bracket_atom -> '[' metal_symbol ']'
bracket_atom -> '[' BAI ']'
bracket_atom -> '[' BAI ringbond ']'

# Defining allowed metal symbols
metal_symbol -> 'Cd' | 'Os' | 'Ti' | 'Rh' | 'Ce' | 'Hg' | 'Cf' | 'Pt' | 'Au' | 'Lu' | 'Cm' | 'Ni' | 'Ho' | 'Nd' | 'Np' | 'Pu' | 'Yb' | 'Tb' | 'Pa' | 'Ag' | 'V' | 'La' | 'U' | 'Ru' | 'Eu' | 'Pd' | 'Zn' | 'Cr' | 'Sm' | 'Am' | 'Dy' | 'Nb' | 'Re' | 'W' | 'Th' | 'Pr' | 'Hf' | 'Sc' | 'Tm' | 'Mn' | 'Mo' | 'Y' | 'Er' | 'Co' | 'Tc' | 'Gd' | 'Zr' | 'Ta' | 'Fe' | 'Bk' | 'Cu' | 'Ir'

# Existing rules for organic and inorganic atoms
aliphatic_organic -> 'B' | 'C' | 'F' | 'H' | 'I' | 'N' | 'O' | 'P' | 'S' | 'Cl' | 'Br' | 'Si' | 'Se'

aromatic_organic -> 'b' | 'c' | 'n' | 'o' | 'p' | 's' | 'se'

# Original BAI and BAC definitions for atoms with isotopes or chirality
BAI -> isotope symbol BAC
BAI -> symbol BAC
BAI -> isotope symbol
BAI -> symbol
BAC -> chiral BAH
BAC -> BAH
BAC -> chiral
BAH -> hcount BACH
BAH -> BACH
BAH -> hcount
BACH -> charge

symbol -> aliphatic_organic
symbol -> aromatic_organic

# Isotopes and valency rules
isotope -> DIGIT
isotope -> DIGIT DIGIT
isotope -> DIGIT DIGIT DIGIT
DIGIT -> '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' | '0'

chiral -> '@' | '@@'
hcount -> 'H' | 'H' DIGIT
charge -> '-' | '-' DIGIT | '-' DIGIT DIGIT | '+' | '+' DIGIT | '+' DIGIT DIGIT

bond -> '-' | '=' | '#' | '/' | '\\'

ringbond -> DIGIT | bond DIGIT

RB -> RB ringbond | ringbond
BB -> BB branch | branch
branch -> '(' chain ')'
branch -> '(' bond chain ')'

Nothing -> None

"""

# form the CFG and get the start symbol
GCFG = nltk.CFG.fromstring(gram)
