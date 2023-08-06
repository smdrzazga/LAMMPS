from math import sqrt
import numpy as np
import functools
from scipy.sparse.linalg import eigsh
import scipy.sparse.linalg

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
        for lines in f:
            i += 1
            if i == 4 :
                return int(lines.strip().split()[0])
            else:
                continue




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
        self.id = int(id)
        self.type = type
        self.position = np.array([float(x), float(y), float(z)])

    def __repr__(self) -> str:
        return f"{self.id} | {self.type} | {self.position}"

    def is_center(self) -> bool:
        _center_id = 6
        _atoms_in_molecule = 11
        if int(self.id) % _atoms_in_molecule == _center_id:
            return True
        else:
            return False
    
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
    def __init__(self, id, num_atoms):
        self.id = id   
        self.comp = []
        self.atoms = num_atoms

    def center(self):
        return self.comp[self.atoms//2]


    def center_of_mass(self):
        if len(self.comp) != self.atoms:
            raise Exception("Error: molecule not fully read." )

        com = np.zeros(3)

        for i in range(self.atoms):
            com += self.comp[i].position
            
        return com / self.atoms


    def polarization(self):
        return self.comp[self.atoms//2].position - (self.comp[-1].position - self.comp[0].position) / 2


    def shift(self, x, y, z):
        m = self.atoms // 2
        mid = Vector(self.comp[m].position[0], self.comp[m].position[1], self.comp[m].position[2] )
        for i in range(self.atoms):     
            self.comp[i].position[0] -= mid.x - x
            self.comp[i].position[1] -= mid.y - y
            self.comp[i].position[2] -= mid.z - z


    def rotate_x(self, theta):
        # change degrees to radians
        theta *= ( np.pi / 180 )
        for i in range(self.atoms):
            temp_y = self.comp[i].position[1]
            temp_z = self.comp[i].position[2]
            
            # self.comp[i].x = self.comp[i].x 
            self.comp[i].position[1] = temp_y * np.cos(theta) - temp_z * np.sin(theta)
            self.comp[i].position[2] = temp_y * np.sin(theta) + temp_z * np.cos(theta)


    def rotate_y(self, theta):
        # change degrees to radians
        theta *= ( np.pi / 180 )
        for i in range(self.atoms):
            temp_x = self.comp[i].position[0]
            temp_z = self.comp[i].position[2]

            self.comp[i].position[0] = temp_x * np.cos(theta) + temp_z * np.sin(theta)
            # self.comp[i].y = self.comp[i].y 
            self.comp[i].position[2] = -1 * temp_x * np.sin(theta) + temp_z * np.cos(theta)

    def rotate_z(self, phi):        
        # change degrees to radians
        phi *= ( np.pi / 180 )
        for i in range(self.atoms):
            temp_x = self.comp[i].position[0]
            temp_y = self.comp[i].position[1]
            
            self.comp[i].position[0] = temp_x * np.cos(phi) - temp_y * np.sin(phi)
            self.comp[i].position[1] = temp_x * np.sin(phi) + temp_y * np.cos(phi)
            # self.comp[i].z = self.comp[i].z 


    def print(self):
        for i in range(Molecule.atoms):
            print(f"{self.comp[i].position}")
        

    def director(self):
        director = self.comp[-1].position - self.comp[0].position
        director /= np.linalg.norm(director)
        
        return director


    def add(self, atom):
        self.comp.append(atom)

class Rescale:
    def __init__(self, current_packing_fraction, target_packing_fraction):
        self.current_packing_fraction = current_packing_fraction
        self.target_packing_fraction = target_packing_fraction

    target_box = Vector(1, 1, 1)
    target_volume = 1
    current_volume = 1
    factor = 1


class Line:
    def __init__(self, id, type, x, y, z):
        self.id = int(id)
        self.type = type
        self.x = float(x)
        self.y = float(y)
        self.z = float(z) 
    
    def is_center(self) -> bool:
        _center_id = 6
        _atoms_in_molecule = 11
        if int(self.id) % _atoms_in_molecule == _center_id:
            return True
        else:
            return False


class Pixel:
    def __init__(self):
        self.components = []

    def assign(self, atom: Atom):
        self.components.append(atom.position)

    def mean(self) -> float:
        # mean = functools.reduce(lambda x,y: x+y, self.components) / len(self.components)
        
        mean = 0
        # check for zero-division
        if len(self.components) > 0:
            mean = np.mean(self.components, axis=0)
            mean /= len(self.components)

        return mean
    
    def colour(self) -> float:
        return 1

    def merge_pixels(self, new_pixel) -> None:
        self.components.extend(new_pixel.components)


class CenterPixel(Pixel):
    def colour(self):
        return len(self.components)


class DirectorPixel(Pixel):
    def __init__(self):
        super().__init__()

    def Q(self):
        Q = np.zeros((3,3))
        
        for director in self.components:
            Q += 0.5 * (3*np.tensordot(director, director, axes=0) - np.identity(3))
        
        if len(self.components) > 0:
            Q /= len(self.components)

        return Q

    def local_director(self):
        if len(self.components) == 0:
            return np.zeros(3)
        
        eigenvalue, eigenvector = scipy.sparse.linalg.eigsh(self.Q(), k=1, which="LM")
        eigenvector = np.reshape(eigenvector, -1)

        if eigenvector[-1] < 0:
            eigenvector *= -1

        return eigenvector


    def P2Value(self, direc=[1,1,1], trunc=True):
        direction = np.array(direc) / np.linalg.norm(np.array(direc))

        # if pixel is empty then return any value
        if len(self.components) == 0:
            return 0
        
        # Legendre 2nd polynomial with respect to dot product (cosine) between average pixel vector and provided direction
        p2value = 0.5 * (3 * ( np.dot(self.local_director(), direction) )**2 - 1)

        if trunc:
            # rescale P2Value from -0.5 to 1, to range 0 to 1
            p2value = (p2value + 0.5) / 1.5 

        return p2value
    
    def colour(self):
        # return self.P2Value()
        return self.local_director()[1]
    

class Screen:
    def __init__(self, x, y, pixel_class):
        self.x = x
        self.y = y
        self.screen = [[pixel_class() for j in range(x)] for i in range(y)]

    def __repr__(self) -> str:
        rows = str()
        for i in range(self.y):
            for j in range(self.x):
                cell = f"{len(self.screen[i][j].components)} "
                rows = rows + cell
            rows = rows + "\n"
        return rows

    def determine_pixel(self, atom: Atom, box: Simulation_box, plane='xz') -> tuple[int, int]:
        # prepare data from proper directions for further processing 
        match plane:
            case "xz":
                box1 = box.x
                box2 = box.z
                atom_position1 = atom.position[0]
                atom_position2 = atom.position[2]

            case "xy":
                box1 = box.x; box2 = box.y
                atom_position1 = atom.position[0]
                atom_position2 = atom.position[1]
        
            case "yz":
                box1 = box.y
                box2 = box.z
                atom_position1 = atom.position[1]
                atom_position2 = atom.position[2]
    
        # determine size of single pixel in the same units as box' length
        self._bin_size_1 = box1 / self.x
        self._bin_size_2 = box2 / self.y
        
        # compute which pixel corresponds to molecules position
        bin1 = int(atom_position1 // self._bin_size_1)
        bin2 = int(atom_position2 // self._bin_size_2)

        # handle exceptions (molecule center does not have to be inside box even when middle atom is inside)
        if bin1 >= self.x: bin1 = self.x - 1
        if bin2 >= self.y: bin2 = self.y - 1
        if bin1 < 0: bin1 = 0
        if bin2 < 0: bin2 = 0

        return bin1, bin2

    def assign(self, atom: Atom, x, y):
        self.screen[y][x].assign(atom)

    def append_screenshot(self, screenshot: list):
        for i in range(self.y):
            for j in range(self.x):
                self.screen[i][j].merge_pixels(screenshot.screen[i][j])
    
    def colour(self):
        colour = np.zeros((self.x, self.y))

        for i in range(self.y):
            for j in range(self.x):
                colour[i, j] = self.screen[i][j].colour()

        return colour

    def avg_director(self):
        if not isinstance(self.screen[0][0], DirectorPixel):
            raise Exception("Pixel has to be Director Pixel")

        Q = np.zeros((3,3))
        for i in range(self.y):
            for j in range(self.x):
                Q += self.screen[i][j].Q()

        Q /= (self.x * self.y)**2
        _, avg_director = scipy.sparse.linalg.eigsh(Q, k=1, which="LA")
        avg_director = np.reshape(avg_director, -1)

        return avg_director


class Screenshot(Screen):
    def __init__(self, x, y, pixel_class):
        super().__init__(x, y, pixel_class)


    def pixels_to_scroll(self, pixel_num, box: Simulation_box, drift) -> int:
        # determine numer of pixels that the image should be scrolled in order to get smectic interferences at the same place
        return int(pixel_num * drift / box.z)
    
    
    def scroll(self, pixels_to_scroll: int) -> None:
        # scrolling as slicing screenshot into two pieces and changing their order 
        self.screen = self.screen[pixels_to_scroll:] + self.screen[:pixels_to_scroll]


class HorizontalSlice(Screen):
    # initializing Slice as a single pixel of type matching current screen
    def __init__(self, screen: Screen, row):
        self.component = type(screen.screen[0][0])()
        self.row = row

    def read_slice(self, screen, start, end):
        for i in range(start, end):
            self.component.merge_pixels(screen.screen[self.row][i])
        
    def slice_director(self):
        if not isinstance(self.slice, DirectorPixel):
            raise Exception("Pixel has to be Director Pixel")

        return self.component.local_director()


class VerticalSlice(Screen):
    # initializing Slice as a single pixel of type matching current screen
    def __init__(self, screen: Screen, row):
        self.component = type(screen.screen[0][0])()
        self.row = row

    def read_slice(self, screen, start, end):
        for i in range(start, end):
            self.component.merge_pixels(screen.screen[i][self.row])
        
    def slice_director(self):
        if not isinstance(self.slice, DirectorPixel):
            raise Exception("Pixel has to be Director Pixel")

        return self.component.local_director()


class SmecticParameter():
    def __init__(self, periods) -> None:
        self.parameter = 0 + 0j
        self.count = 0
        self.periods = 4

    def __repr__(self) -> str:
       return f"{np.abs(self.parameter):.3f} | {self.count}"

    # only for vertical slices in y-z plane!
    def add_atom(self, center, box: Simulation_box):
        self.parameter += np.exp(self.periods * 2*np.pi*1j * center[-1] / box.z)
        self.count += 1
        # print(np.exp(2*np.pi*1j * center[-1] / box.z))

    def normalize(self):
        if self.count != 0:
            self.parameter /= self.count
    
    def read_screen(self, screen: Screen, start, end, box: Simulation_box):
        for x in range(start, end):
            for y in range(screen.y):
                # for coords in screen.screen[y][x].components:
                #     self.add_atom(coords, box)
                self.parameter += np.exp(self.periods * 2*np.pi*1j * y / screen.y) * screen.screen[y][x].colour()
                self.count += screen.screen[y][x].colour()

def wrap_atom_to_box(atom: Atom, box: Simulation_box) -> Atom:
    atom.position[0] = atom.position[0] % box.x
    atom.position[1] = atom.position[1] % box.y
    atom.position[2] = atom.position[2] % box.z

    return atom

        

