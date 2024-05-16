from __future__ import print_function, absolute_import
import os, sys
from rdkit import Chem as C
from rdkit.Chem import Descriptors as D
from rdkit.Chem import rdMolDescriptors as CD

def get_stuff(smiles):

        mol = C.MolFromSmiles(smiles)
        #hac = D.HeavyAtomCount(mol)

        mw = CD.CalcExactMolWt(mol)
        logp = C.Crippen.MolLogP(mol)
        rotB = D.NumRotatableBonds(mol)
        HBA = CD.CalcNumHBA(mol)
        HBD = CD.CalcNumHBD(mol)
        q = C.GetFormalCharge(mol)

        print(mw, logp, rotB, HBA, HBD)
        print("MW is ",mw)
        print("logP is ",logp)
        print("Rotatable Bonds is ",rotB)
        print("HB Donors is ",HBD)
        print("HB Acceptors is ",HBA)
        print("Formal Charge is ",q)
        return(mw, logp, rotB, HBA, HBD, q)


def main():

        pwd = os.getcwd()+"/"
        smiles = sys.argv[1]

        get_stuff(smiles)




main()
