import sys
from openmoltools import amber_parser

if __name__ == "__main__":
    if len(sys.argv[1:]) <= 1:
        print("""Usage: processAmberForceField.py some_path/gaff.dat ligand_name.mol2 ligand_name.frcmod
""")
    else:
        parser = amber_parser.AmberParser()
        parser.parse_filenames(sys.argv[1:])
        stream = parser.generate_xml()
        print(stream.read())

