
# Code for reference data:
                
pt = Chem.GetPeriodicTable()
atoms = ['H', 'He',
         'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
         'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
         'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
         'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe']
covRadii = {}
masses = {}
for atom in atoms:
    covRadii[atom] = pt.GetRcovalent(atom)
    masses[atom] = pt.GetAtomicWeight(atom)
print(covRadii)
print(masses) # Then copy pasted both dictionaries into notebook

from rdkit import Chem
from rdkit.Chem import AllChem


# Code which was used to generate .xyz files for this project:

m1 = Chem.MolFromSmiles('CC') # Replace SMILES in this line for desired molecule
m1 = Chem.AddHs(m1)
AllChem.EmbedMolecule(m1) # Generate 3D coords
AllChem.MMFFOptimizeMolecule(m1)
def write_xyz_file(fragment, fragment_name):
    ''' Adapted from: https://github.com/jensengroup/take_elementary_step/blob/master/write_input_files.py '''
    number_of_atoms = fragment.GetNumAtoms()
    symbols = [a.GetSymbol() for a in fragment.GetAtoms()] 
    for idx,conf in enumerate(fragment.GetConformers()):
        file_name = fragment_name+".xyz"
        with open(file_name, "w") as file:
            file.write(str(number_of_atoms)+"\n")
            file.write("\n")
            for atom,symbol in enumerate(symbols):
                p = conf.GetAtomPosition(atom)
                line = " ".join((symbol,str(p.x),str(p.y),str(p.z),"\n"))
                file.write(line)
write_xyz_file(m1,'CH3CH3')



