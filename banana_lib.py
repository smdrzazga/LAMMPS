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
    
    
    def shift(self, x, y, z):
        m = self.atoms // 2
        mid = Vector(self.comp[m].x, self.comp[m].y, self.comp[m].z )
        for i in range(self.atoms):     
            self.comp[i].x -= mid.x - x
            self.comp[i].y -= mid.y - y
            self.comp[i].z -= mid.z - z


    def rotate_x(self, theta):
        # change degrees to radians
        theta *= ( np.pi / 180 )
        for i in range(self.atoms):
            temp_y = self.comp[i].y
            temp_z = self.comp[i].z
            
            # self.comp[i].x = self.comp[i].x 
            self.comp[i].y = temp_y * np.cos(theta) - temp_z * np.sin(theta)
            self.comp[i].z = temp_y * np.sin(theta) + temp_z * np.cos(theta)


    def rotate_y(self, theta):
        # change degrees to radians
        theta *= ( np.pi / 180 )
        for i in range(self.atoms):
            temp_x = self.comp[i].x
            temp_z = self.comp[i].z

            self.comp[i].x = temp_x * np.cos(theta) + temp_z * np.sin(theta)
            # self.comp[i].y = self.comp[i].y 
            self.comp[i].z = -1 * temp_x * np.sin(theta) + temp_z * np.cos(theta)

    def rotate_z(self, phi):        
        # change degrees to radians
        phi *= ( np.pi / 180 )
        for i in range(self.atoms):
            temp_x = self.comp[i].x
            temp_y = self.comp[i].y
            
            self.comp[i].x = temp_x * np.cos(phi) - temp_y * np.sin(phi)
            self.comp[i].y = temp_x * np.sin(phi) + temp_y * np.cos(phi)
            # self.comp[i].z = self.comp[i].z 

    def print_molecule(self):
        for i in range(Molecule.atoms):
            print("%.5f, %.5f, %.5f" % (self.comp[i].x, self.comp[i].y, self.comp[i].z))

class Rescale:
    def __init__(self, current_packing_fraction, target_packing_fraction):
        self.current_packing_fraction = current_packing_fraction
        self.target_packing_fraction = target_packing_fraction

    target_box = Vector(1, 1, 1)
    target_volume = 1
    current_volume = 1
    factor = 1

