from math import sqrt
import numpy as np

def is_molecule_torn(molecule, axis):
    
    is_torn = 0
    if(axis == 'x'):
        for i in range(len(molecule)-1):
            if abs(molecule[i].x - molecule[i+1].x ) > 10:
                is_torn = 1    


    if(axis == 'y'):
        for i in range(len(molecule)-1):
            if abs(molecule[i].y - molecule[i+1].y ) > 10:
                is_torn = 1    

    if(axis == 'z'):
        for i in range(len(molecule)-1):
            if abs(molecule[i].z - molecule[i+1].z ) > 10:
                is_torn = 1    

    return is_torn


def sew_molecule(molecule, x_box, y_box, z_box):
    # check whether whole 11 atom molecule is read
    if len(molecule.comp) != 11: 
        return 

    # sewing back torn atoms by moving them box size in positive direction  
    if is_molecule_torn(molecule.comp, "x") == 1: 
        for j in range(11):
            if molecule.comp[j].x < 10:
                molecule.comp[j].x += x_box
    if is_molecule_torn(molecule.comp, "y") == 1: 
        for j in range(11):    
            if molecule.comp[j].y < 10:
                molecule.comp[j].y += y_box
    if is_molecule_torn(molecule.comp, "z") == 1: 
        for j in range(11):                
            if molecule.comp[j].z < 10:
                molecule.comp[j].z += z_box    


def write_heading(target, n, trim):
    
    # writing molecule file heading
    with open (target, "w") as t:
        t.write("LAMMPS Description\n \n") 

        t.write("%d atoms\n" % n)
        t.write("0 bonds\n")
        t.write("0 angles\n")
        t.write("0 dihedrals\n")
        t.write("0 impropers\n\n")

        t.write("1 atom types\n\n")

        print("0 %f xlo xhi" % trim.x, file=t)
        print("0 %f ylo yhi" % trim.y, file=t)
        print("0 %f zlo zhi \n" % trim.z, file=t)

        t.write("Masses\n\n")

        t.write("1 1\n\n")

        t.write("Atoms \n\n")



def read_boundaries(file):

    i = 0    

    with open (file, "r") as f:
        for line in f:
            
            i += 1
            
            if i == 6:
                x_min = float(line.strip().split()[0])
                x_max = float(line.strip().split()[1])
                continue
            if i == 7:
                y_min = float(line.strip().split()[0])
                y_max = float(line.strip().split()[1])
                continue    
            if i == 8:
                z_min = float(line.strip().split()[0])
                z_max = float(line.strip().split()[1])
                continue    

            if i > 10: 
                return x_min, x_max, y_min, y_max, z_min, z_max


def read_number_of_atoms(location):
    with open(location, "r") as f:

        i = 0
        for line in f:
            l = line.strip().split()
            i += 1

    return int(l[0])


def find_max(file, coord):
    
    max = 0.
    rows = {"x": 3, "y": 4, "z": 5}
    i = 0

    with open (file, "r") as f:
        for line in f:

            i += 1

            # skip headline
            if i < 30: 
                continue
        
            # compare x coordinate
            l = line.strip().split()
            x = float(l[rows[coord]])

            if (max < x):
                max = x

    return max


def find_min(file, coord):
    
    min = 1000.
    rows = {"x": 3, "y": 4, "z": 5}
    i = 0

    with open (file, "r") as f:
        for line in f:

            i += 1

            # skip headline
            if i < 30: 
                continue
        
            # compare x coordinate
            l = line.strip().split()
            x = float(l[rows[coord]])

            if (x < min):
                min = x

    return min


def scale(target_volume, current_volume, num_walls):
    # calculate ratio between current and new box edges lengths
    if num_walls < 3:
        return pow(target_volume / current_volume, 1 / (3 - num_walls))
    if num_walls == 3:    
        return 1


class Atom:
    def __init__(self, id, type, x, y, z):
        self.id = id
        self.type = type
        self.x = x
        self.y = y
        self.z = z
    
    volume = 1.333 * 3.14 * 0.5 * 0.5 * 0.5


class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


class Simulation_box:
    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max, atoms):
        self.min = Vector(x_min, y_min, z_min)
        self.max = Vector(x_max, y_max, z_max)
        self.x = x_max - x_min
        self.y = y_max - y_min
        self.z = z_max - z_min

        self.atoms = atoms
        self.mols = 1

        self.volume = self.x * self.y * self.z

class Molecule:
    def __init__(self, id):
        self.id = id   
        self.comp = []

    atoms = 11
    def center_of_mass(self):
        if len(self.comp) != self.atoms:
            return "Error: molecule not fully read." 

        x_com = 0
        y_com = 0
        z_com = 0

        for i in range(self.atoms):
            x_com += self.comp[i].x
            y_com += self.comp[i].y
            z_com += self.comp[i].z
            
        return Vector(x_com/self.atoms, y_com/self.atoms, z_com/self.atoms)



class Rescale:
    def __init__(self, current_packing_fraction, target_packing_fraction):
        self.current_packing_fraction = current_packing_fraction
        self.target_packing_fraction = target_packing_fraction

    target_box = Vector(1, 1, 1)
    target_volume = 1
    current_volume = 1
    factor = 1



# location of snapshot file and path to output file
location = "C:/Users/komok/Desktop/final_snapshot_0.32.lammpstrj"
# location = "C:/Users/komok/Desktop/tworzenie_atomow_z_pliku.txt"
target = "C:/Users/komok/Desktop/tworzenie_atomow_z_pliku.txt"

target_packing_fraction = 0.31

# read box size
box = Simulation_box( *read_boundaries(location), read_number_of_atoms(location))
box.mols = box.atoms / Molecule.atoms

# read packing fraction and the one to be obtained after rescaling
current_packing_fraction = box.atoms * Atom.volume / box.volume 
rescale = Rescale(current_packing_fraction, target_packing_fraction)



# target box volume
rescale.target_volume = box.atoms * Atom.volume / rescale.target_packing_fraction
rescale.current_volume = box.volume
rescale.factor = scale(rescale.target_volume, rescale.current_volume, 1)
# calculate scale factor, by which sides y-z are to be multiplied in order to get target packing fraction
rescale.target_box = Vector(box.x, box.y * rescale.factor, box.z * rescale.factor)



# creating middle step file with molecules wrapped back to original box and sewed if torn
with open (location, "r") as f:
    with open (target, "w") as t:

        i = 0
        for line in f:
            i += 1
            elements = line.strip().split()

            if len(elements) != 5:
                continue
            try: 
                test = int(elements[0])   
            except ValueError:
                continue
            
            atom = Atom( int(elements[0]), elements[1], float(elements[2]), float(elements[3]), float(elements[4]) )

            # atom id and coordinates
            atom.id = int(elements[0])   
            atom.type = elements[1]
            atom.x = float(elements[2])
            atom.y = float(elements[3])
            atom.z = float(elements[4])  

            # which molecule atom belongs to
            molecule_id = int(1 + (atom.id - 1)/Molecule.atoms)

            # box origin shifted to (0,0,0) and atom coordinates wrapped to be in that box
            atom.x = (atom.x - box.min.x) % (box.x)
            atom.y = (atom.y - box.min.y) % (box.y)
            atom.z = (atom.z - box.min.z) % (box.z)

            # for every molecule its atom coordinates are temporary stored in a list, cleared every 11 atoms starting from the first one
            if atom.id % 11 == 1:
                molecule = Molecule(molecule_id)

            molecule.comp.append(atom)

            # sewing back torn atoms by moving them box size in positive direction  
            sew_molecule(molecule, box.x, box.y, box.z)

            # writing wrapped, sewed atom coordinates to target file 
            if len(molecule.comp) == 11:    

                for j in range(11):
                
                    t.write('{} '.format(molecule.comp[j].id))
                    t.write('{} '.format(molecule.id))
                    t.write('{} '.format(molecule.comp[j].type))
                    t.write('{} '.format(molecule.comp[j].x + 1.5))
                    t.write('{} '.format(molecule.comp[j].y))
                    t.write('{} \n'.format(molecule.comp[j].z))




# finding largest x,y,z coordinates arcross all the wrapped, sewed atoms
xyz_max = Vector(find_max(target, "x"), find_max(target, "y"), find_max(target, "z"))

write_heading(target, box.atoms, rescale.target_box)
with open (location, "r") as f:
    with open (target, "a") as t:
    
        i = 0
        
        for line in f:
            # if i < 10: 
            #     continue
            i += 1
            elements = line.strip().split()
            
            if len(elements) != 5:
                continue

            try: 
                test = int(elements[0])   
            except ValueError:
                continue


            # atom id and coordinates
            atom = Atom( int(elements[0]), elements[1], float(elements[2]), float(elements[3]), float(elements[4]) )

            # which molecule atom belongs to
            molecule_id = int(1 + (atom.id - 1)/Molecule.atoms)

            # box origin shifted to (0,0,0) and atom coordinates wrapped to be in that box
            atom.x = (atom.x - box.min.x) % (box.x)
            atom.y = (atom.y - box.min.y) % (box.y)
            atom.z = (atom.z - box.min.z) % (box.z)

            # for every molecule its atom coordinates are temporary stored in a list, cleared every 11 atoms starting from the first one
            if atom.id % 11 == 1:
                molecule = Molecule(molecule_id)

            molecule.comp.append(atom)

            # sewing back torn atoms by moving them box size in positive direction  
            sew_molecule(molecule, box.x, box.y, box.z)

            # writing wrapped, sewed atom coordinates to target file 
            if len(molecule.comp) == Molecule.atoms:    

                for j in range(11):
                
                    print(molecule.comp[j].id, end = " ", file=t)
                    print(molecule.id, end = " ", file=t)
                    print(1, end = " ", file=t)
                    print("%f " % (molecule.comp[j].x - molecule.center_of_mass().x * (1 - (rescale.target_box.x - 3) / xyz_max.x ) + 1.5), end=" ", file=t)
                    print("%f " % (molecule.comp[j].y - molecule.center_of_mass().y * (1 - (rescale.target_box.y - 3.5) / xyz_max.y ) + 1.5), end=" ", file=t)
                    print("%f " % (molecule.comp[j].z - molecule.center_of_mass().z * (1 - (rescale.target_box.z - 3.5) / xyz_max.z ) + 1.5), file=t)



print("\nMaximum x,y,z coordinates after sewing molecules:")
print(find_max(target, "x"))
print(find_max(target, "y"))
print(find_max(target, "z"))

print("\nMinimum x,y,z coordinates after sewing molecules:")
print(find_min(target, "x"))
print(find_min(target, "y"))
print(find_min(target, "z"))

# box size
print("\nOriginal box size:")
print(box.x)
print(box.y)
print(box.z)
print("\nTarget box size:")
print(rescale.target_box.x)
print(rescale.target_box.y)
print(rescale.target_box.z)

print("\nOriginal volume: ", box.x * box.y * box.z)
print("Target volume: ", rescale.target_volume)
print("Scale: ", rescale.factor)

print("\nNumber of atoms:")
print(box.atoms)

print("\nNumber of molecules:")
print(box.mols)

print("\nOriginal packing fraction: %f" % rescale.current_packing_fraction)
print("Target packing fraction: %f\n" % rescale.target_packing_fraction)

