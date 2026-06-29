import math
from cov_radii import cov_dictionary

# Function: parse_line
# This function reads a line from a molecular coordinate file and extracts atom information.
def parse_line(number, string):
    parts = string.split() 
    element_symbol = parts[0]  
    x, y, z = map(float, parts[1:])  
    return [number, element_symbol, [x, y, z]]

# Function: distance
# This function calculates the Euclidean distance between two 3D points.
def distance(list1, list2):
    try:
        return math.sqrt(
            (list2[0]- list1[0])**2 +
            (list2[1]- list1[1])**2 +
            (list2[2]- list1[2])**2
        )
    except (IndexError, TypeError, NameError) as e:
        print(f"Error: {e}")
        return None

# Function: coordination
# This function determines which atoms are bonded to a given atom.
# It compares the distance between the given atom (list1) and every other atom in list2.
def coordination(list1, list2):
    bonded_atoms = []
    for atom in list2:
        if atom == list1:
            continue
        r1 = cov_dictionary[list1[1]]
        r2 = cov_dictionary[atom[1]]
        d = distance(list1[2], atom[2])
        if d <= r1 + r2 + 0.4:
            bonded_atoms.append(atom)
    return bonded_atoms

# Function: is_acid
# This function checks if a given carbon atom is part of a carboxylic acid group.
# It only applies if the central atom (list1) is a carbon ("C").
# The function looks at all bonded atoms (list2). If at least two oxygen ("O") atoms
# are directly bonded to the carbon, it classifies the atom as being in an acid group.
def is_acid(list1, list2):
    oxygens = []
    if list1[1] != "C": 
        return False
    for atom in list2:
        if atom[1] == "O":
            oxygens.append(atom)
    if len(oxygens) >= 2:
        oxygens.append(list1)
        return True
    else:
        return False

# Function: is_terminal
# This function checks whether a carbon atom is a terminal carbon.
# A terminal carbon is defined here as a carbon atom bonded to either:
# - only one other carbon atom, OR
# - no other carbon atoms at all.
# Returns True if the carbon is terminal, False otherwise.
def is_terminal(list1, list2):
    if list1[1] != "C":
        return False 
    carbon_count = 0
    for atom in list2:
        if atom[1] == "C":
            carbon_count += 1
    if carbon_count == 1:
        return True
    elif carbon_count == 0:
        return True
    else:
        return False

# Function: unsaturation
# This function estimates the degree of unsaturation for a carbon atom.
# It looks at the number and type of atoms bonded to the carbon.
def unsaturation(list1, list2):
    if list1[1] != 'C':
        return None

    saturation_bucket = [atom[1] for atom in list2]

    if len(saturation_bucket) == 4:
        return 0  

    elif len(saturation_bucket) == 3 and saturation_bucket.count('C') == 1 and saturation_bucket.count('O') == 2:
        return 1  

    elif len(saturation_bucket) == 3:
        return 1  

    elif len(saturation_bucket) == 2:
        return 2  

    else:
        return None

# Function: smiles
# This function generates a simplified SMILES (Simplified Molecular Input Line Entry System) string
# representation of the molecule.
def smiles(atoms):
    smiles_str = ""

    start_atom = None
    start_neighbours = None
    for atom in atoms:
        if atom[1] == "C":
            bonds = coordination(atom, atoms)
            if is_acid(atom, bonds):
                smiles_str += "OC(=O)"
                start_atom = atom
                start_neighbours = bonds
                break
    
    if start_atom is None:
        return smiles_str  
    
    next_atom = None
    for n in start_neighbours:
        if n[1] == "C":
            next_atom = n
            next_neighbours = coordination(n, atoms)
            break
    
    if next_atom is None:
        return smiles_str  

    prev_atom = start_atom
    prev_neighbours = start_neighbours
    curr_atom = next_atom
    curr_neighbours = next_neighbours

    while True:
        u_prev = unsaturation(prev_atom, prev_neighbours)
        u_curr = unsaturation(curr_atom, curr_neighbours)

        if u_curr == 0:
            smiles_str += "C"
        elif u_prev == 1 and u_curr == 1:
            smiles_str += "C"
        elif u_prev == 0 and u_curr == 1:
            smiles_str += "C="
        elif u_prev == 0 and u_curr == 2:
            smiles_str += "C#"
        else:
            smiles_str += "C"

        if is_terminal(curr_atom, curr_neighbours):
            break

        track_prev = prev_atom
        prev_atom = curr_atom
        prev_neighbours = curr_neighbours

        next_atom = None
        for n in curr_neighbours:
            if n[1] == "C" and n != track_prev:
                next_atom = n
                next_neighbours = coordination(n, atoms)
                break

        if next_atom is None:
            smiles_str += "C"
            break

        curr_atom = next_atom
        curr_neighbours = next_neighbours

    return smiles_str

# Function: main
# It combines all previous functions to generate a SMILES string from a given input file
# The file must be an XYZ file format
def main():
    file_path = input("Enter a file: ")
    atoms = []
    try:
        with open(file_path, "r") as file:
            lines = file.readlines()
        i = 1  
        for line in lines[2:]: 
            atom = parse_line(i, line)
            atoms.append(atom)
            i += 1
        smiles_string = smiles(atoms)
        print("SMILES:", smiles_string)

    except FileNotFoundError:
        print("This file was not found")
if __name__ == "__main__":
    main()