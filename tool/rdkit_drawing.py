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


def mol2svg(mol,molSize=(400,400),kekulize=True, **kwargs):
    '''
    Create svg figure of a mol

    Parameters
    ----------
    mol : mol
        mol to draw

    molSize : tuple
        size of the figure

    kekulize : bool
        see rdkit docu

    Returns
    -------
    
    SVG : string
    Figure of the svg

    '''

    # mc = Chem.Mol(mol.ToBinary())
    # print(mc)

    # if kekulize:
    #     try:
    #         Chem.Kekulize(mc)
    #     except:
    #         mc = Chem.Mol(mol.ToBinary())
    # if not mc.GetNumConformers():
    #     rdDepictor.Compute2DCoords(mc)

    # print(Chem.MolToSmiles(mol))

    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mol, **kwargs)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg.replace('svg:','')

def smiles2svg(smiles, **kwargs):

    mol = Chem.MolFromSmiles(smiles)

    svg = mol2svg(mol, **kwargs)

    return(svg)