####################
#Rdkit plotting
####################

from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import DataStructs
import os
import numpy as np
import collections

import itertools
# from IPython.display import display, SVG

def DrawBitFromSmiles(smiles, fp_choice):

    fp_choice = int(fp_choice)

    mol = Chem.MolFromSmiles(smiles)

    bitInfo = {}
    fp = AllChem.GetMorganFingerprint(mol, 
                                      3, 
                                      bitInfo=bitInfo, 
                                      useChirality=True)
    
    svg = Draw.DrawMorganBit(mol, fp_choice, bitInfo, useSVG=True)
    svg = svg.replace("\n","")
    # exit()
    return(svg)