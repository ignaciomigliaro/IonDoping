# Ion Doping
This repo serves as a way to randomly dope Crystal Structures with anions and balance charge with removal of cations. This is to account that the functions provided by ASE or Pymatgen are desgined for doping cations. This is designed to run for VASP, and for now the cation balance charge assumes the ion is of charge -1. 

usage: doping.py [-h] [--target_atom TARGET_ATOM] [--anion ANION] [--cation CATION] [--fraction FRACTION] [--count COUNT] input_files output_dir

This is a code to generate doped structures from an undoped input file (cif,poscar)

positional arguments:
  input_files           Input files of the VASP file containing POSCAR of undoped structure, KPOINTS, POTCAR,INCAR
  output_dir            Directory where the input files for either QE or VASP will be written

options:
  -h, --help            show this help message and exit
  --target_atom TARGET_ATOM
                        Element to be replaced
  --anion ANION         Element which you wish to dope with
  --cation CATION       Cation to balance the electrolyte charge
  --fraction FRACTION   Fraction of target atoms being replaced
  --count COUNT         Number of doped structures you want generated