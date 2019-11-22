import pandas
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import openbabel
import numpy as np 
import os 
from biopandas.mol2 import PandasMol2
import openbabel
import argparse

def SYBYL(row):
    try:
        obconversion = openbabel.OBConversion()
        obconversion.SetInAndOutFormats("smi", "mol2")
        mol = openbabel.OBMol()
        obconversion.ReadString(mol,row["SMILES"])  # read molecule from database 
        mol.AddHydrogens()
        output_mol2 = obconversion.WriteString(mol)  # transform smiles into mol2
        file = open("molecule.mol2","w+")   # write mol2 format into the file, molecule.mol2.
        file.write(output_mol2)
        file.close()
        molecule_mol2 = PandasMol2().read_mol2("molecule.mol2")  # use biopandas to static the discriptors
        for element in molecule_mol2.df['atom_type'].value_counts().index:
            if element == 'Al':
                row['Al'] = molecule_mol2.df['atom_type'].value_counts()['Al']
            if element == 'B':
                row['B'] = molecule_mol2.df['atom_type'].value_counts()['B']
            if element == 'Br':
                row['Br'] = molecule_mol2.df['atom_type'].value_counts()['Br']
            if element == 'C.1':
                row['C.1'] = molecule_mol2.df['atom_type'].value_counts()['C.1']
            if element == 'C.2':
                row['C.2'] = molecule_mol2.df['atom_type'].value_counts()['C.2']
            if element == 'C.3':
                row['C.3'] = molecule_mol2.df['atom_type'].value_counts()['C.3']
            if element == 'C.ar':
                row['C.ar'] = molecule_mol2.df['atom_type'].value_counts()['C.ar']
            if element == 'C.cat':
                row['C.cat'] = molecule_mol2.df['atom_type'].value_counts()['C.cat']
            if element == 'Ca':
                row['Ca'] = molecule_mol2.df['atom_type'].value_counts()['Ca']
            if element == 'Cl':
                row['Cl'] = molecule_mol2.df['atom_type'].value_counts()['Cl']
            if element == 'F':
                row['F'] = molecule_mol2.df['atom_type'].value_counts()['F']
            if element == 'H':
                row['H'] = molecule_mol2.df['atom_type'].value_counts()['H']  
            if element == 'Li':
                row['Li'] = molecule_mol2.df['atom_type'].value_counts()['Li']  
            if element == 'Mg':
                row['Mg'] = molecule_mol2.df['atom_type'].value_counts()['Mg']        
            if element == 'N.1':
                row['N.1'] = molecule_mol2.df['atom_type'].value_counts()['N.1']         
            if element == 'N.2':
                row['N.2'] = molecule_mol2.df['atom_type'].value_counts()['N.2']         
            if element == 'N.3':
                row['N.3'] = molecule_mol2.df['atom_type'].value_counts()['N.3']         
            if element == 'N.4':
                row['N.4'] = molecule_mol2.df['atom_type'].value_counts()['N.4']         
            if element == 'N.am':
                row['N.am'] = molecule_mol2.df['atom_type'].value_counts()['N.am']         
            if element == 'N.ar':
                row['N.ar'] = molecule_mol2.df['atom_type'].value_counts()['N.ar']        
            if element == 'N.pl3':
                row['N.pl3'] = molecule_mol2.df['atom_type'].value_counts()['N.pl3']         
            if element == 'Na':
                row['Na'] = molecule_mol2.df['atom_type'].value_counts()['Na']        
            if element == 'O.2':
                row['O.2'] = molecule_mol2.df['atom_type'].value_counts()['O.2']        
            if element == 'O.3':
                row['O.3'] = molecule_mol2.df['atom_type'].value_counts()['O.3'] 
            if element == 'O.co2':
                row['O.co2'] = molecule_mol2.df['atom_type'].value_counts()['O.co2']  
            if element == 'P.3':
                row['P.3'] = molecule_mol2.df['atom_type'].value_counts()['P.3']    
            if element == 'S.2':
                row['S.2'] = molecule_mol2.df['atom_type'].value_counts()['S.2']        
            if element == 'S.3':
                row['S.3'] = molecule_mol2.df['atom_type'].value_counts()['S.3']
            if element == 'S.O2':
                row['S.O2'] = molecule_mol2.df['atom_type'].value_counts()['S.O2']   
            if element == 'S.O':
                row['S.O'] = molecule_mol2.df['atom_type'].value_counts()['S.O']       
            if element == 'Si':
                row['Si'] = molecule_mol2.df['atom_type'].value_counts()['Si']    
            if element == 'Zn':
                row['Zn'] = molecule_mol2.df['atom_type'].value_counts()['Zn']       
        return row
    except:
        print(row["SMILES"],"SYBYL something is wrong!!")


def ECFP(row):
    try:
        mol = Chem.MolFromSmiles(row['SMILES'])
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits).ToBitString()
        for i in range(nbits):
            row[str(i)] = fp[i]
        return row  
    except:
         print(row["SMILES"],"ECFP feature something is wrong!!")


def smile2ECFP_number(row):
    from rdkit.Chem import AllChem
    from rdkit import Chem
    from rdkit.Chem import Draw
    import pandas
    from collections import Counter
    m = Chem.MolFromSmiles(row['SMILES'])
    
    info={}
    fp = AllChem.GetMorganFingerprintAsBitVect(m,radius,nBits=nbits,bitInfo=info).ToBitString()
    
    try:
        for bit in info:
            row["num_"+str(bit)] = len(info[bit])
        return row  
    except:
         print(row["SMILES"],"ECFP feature something is wrong!!")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--smile_file',
                        help='experiment file', default=None)
    parser.add_argument('-n', '--nbits',
                        help='how many bits of ECFP', default=2048, type=int)
    parser.add_argument('-r', '--radius',
                        help='ECFP radius', default=2, type=int)
    
    
    args = vars(parser.parse_args())


    # load data            
    data = pd.read_csv(args["smile_file"])
    # set parameter of ECFP, namely nbits, radius
    nbits = args["nbits"]
    radius = args["radius"]
    # discriptor ECFP 
    data_ECFP = data.reindex(columns = data.columns.tolist() + [str(i) for i in range(nbits)]) 
    data_ECFP = data_ECFP.apply(ECFP, axis=1)

    data_ECFP.to_csv("ECFP.csv",index=False)

    # discriptor ECFP_number 
    data_ECFP_num = data.reindex(columns = data.columns.tolist() + ["num_"+str(i) for i in range(nbits)]) 
    data_ECFP_num = data_ECFP_num.apply(smile2ECFP_number, axis=1).fillna(0)
    data_ECFP_num.to_csv("ECFP_num.csv",index=False)
    # discriptor SYBYL 
    atom_type=['Al','B','Br','C.1','C.2','C.3','C.ar','C.cat','Ca','Cl','F','H',
        'Li','Mg','N.1','N.2','N.3','N.4','N.am','N.ar','N.pl3','Na',
        'O.2','O.3','O.co2','P.3','S.2','S.3','S.O','S.O2','Si','Zn']

    data_SYBYL = data.reindex(columns = data.columns.tolist() + atom_type) 
    data_SYBYL = data_SYBYL.apply(SYBYL, axis=1).fillna(0)
    data_SYBYL.to_csv("SYBYL.csv",index=False)
    os.remove("./molecule.mol2")
