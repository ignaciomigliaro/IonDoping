import argparse
from ase.io import read,write
import random
import os
import shutil

def parse_args(): 
    parser = argparse.ArgumentParser(
        description="This is a code to generate doped structures from an undoped input file (cif,poscar)"
    )
    parser.add_argument(
        "input_files",
        help="Input files of the VASP file containing POSCAR of undoped structure, KPOINTS, POTCAR,INCAR"
    )
    parser.add_argument(
        "output_dir",
        help='Directory where the input files for either QE or VASP will be written'
    )
    parser.add_argument(
        '--target_atom',
        help='Element to be replaced'
    )
    parser.add_argument(
        '--anion',
        help='Element which you wish to dope with'
    )
    parser.add_argument(
        '--cation',
        help='Cation to balance the electrolyte charge'
    )
    parser.add_argument(
        '--fraction',
        help='Fraction of target atoms being replaced'
    )
    parser.add_argument(
        '--count',
        help="Number of doped structures you want generated"
    )
    return parser.parse_args()

def read_file(file_path):
    try:
        # Check if the file is in POSCAR format
        poscar_structure = Poscar.from_file(file_path).structure
        return poscar_structure
    except Exception as poscar_error:
        # If an error occurs, it might not be in POSCAR format
        print(f"Not a POSCAR file. Error: {poscar_error}")

    try:
        # Check if the file is in CIF format
        cif_parser = CifParser(file_path)
        cif_structure = cif_parser.get_structures()[0]
        return cif_structure
    except Exception as cif_error:
        # If an error occurs, it might not be in CIF format
        print(f"Not a CIF file. Error: {cif_error}")

    # If neither POSCAR nor CIF format is detected
    raise ValueError(f"Unsupported file format for file: {file_path}")

def dope_fraction(structure, old_atom_type, new_atom_type, doping_fraction):
    atoms_to_dope = [atom for atom in structure if atom.symbol == old_atom_type]
    num_atoms_to_dope = int(len(atoms_to_dope) * doping_fraction)

    # Randomly select atoms to dope
    atoms_to_dope = random.sample(atoms_to_dope, num_atoms_to_dope)

    # Replace selected atoms with the new atom type
    for atom in atoms_to_dope:
        atom.symbol = new_atom_type

    return num_atoms_to_dope

def balance_charge(structure, cation_type, num_cations_removed):
    # Remove cations to balance the charge
    removed_cations = 0
    for _ in range(num_cations_removed):
        cation_indices = [i for i, atom in enumerate(structure) if atom.symbol == cation_type]
        if cation_indices:
            random_cation_index = random.choice(cation_indices)
            structure.pop(random_cation_index)
            removed_cations += 1

    return removed_cations

def file_dir_iterate(input_directory,output_directory,anion_target,cation_type,anion_replace,num_structures=10,doping_fraction=.5):
    poscar_file = input_directory + '/POSCAR'
    initial_structure = read(poscar_file)
    common_files = ['KPOINTS', 'POTCAR', 'INCAR']
    os.makedirs(output_dir, exist_ok=True)

    common_files = ['KPOINTS', 'POTCAR', 'INCAR']
    for count in range(1, num_structures + 1):
        structure_directory = os.path.join(output_directory, f'dope{count}')
        os.makedirs(structure_directory, exist_ok=True)

        for common_file in common_files:
            source_path = os.path.join(input_directory, common_file)
            destination_path = os.path.join(structure_directory, common_file)
            shutil.copy(source_path, destination_path)
    
        structure = initial_structure.copy()
        num_replaced_anions = dope_fraction(structure, anion_target, anion_replace, doping_fraction)
        print(f"Structure {count}: Replaced {num_replaced_anions} O atoms with F atoms.")
        num_removed_cations = balance_charge(structure, cation_type, num_replaced_anions)
        print(f"Structure {count}: Removed {num_removed_cations} Li atoms to balance the charge.")

        output_poscar = os.path.join(structure_directory, 'POSCAR')
        write(output_poscar, structure)



args = parse_args()
input_dir = args.input_files
output_dir = args.output_dir
target_anion = str(args.target_atom )
anion_replace = str(args.anion)
cation = str(args.cation)
fraction = float(args.fraction)
num_structures = int(args.count)
fraction = float(args.fraction)

file_dir_iterate(input_dir,output_dir,target_anion,cation,anion_replace,num_structures,fraction)