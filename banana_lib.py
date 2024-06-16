from math import sqrt
import numpy as np
import functools
from scipy.sparse.linalg import eigsh
import scipy.sparse.linalg
from mmap import mmap

def is_molecule_torn(molecule, axis: str) -> bool:
    
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


def sew_molecule(molecule, x_box: float, y_box: float, z_box: float) -> None:
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


def write_heading(target: str, n: int, trim) -> None:
    
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


def get_line(file, is_mmap):
    line = file.readline()
    if is_mmap:
        line = line.decode()
    line = line.split()

    return line

def read_boundaries(file: str|mmap, N_TOTAL, open=True, is_mmap=False) -> tuple:
    if open:
        file = open(file, 'r+')

    i = 0
    line = get_line(file, is_mmap)
                    
    while len(line) < 2 or line[1] != 'BOX':
        i += 1
        line = get_line(file, is_mmap)

        if not line:
            return (0, 0, 0, 0, 0, 0)

        if i > N_TOTAL:
            raise Exception("Box dimensions not found!")
    
    l = file.readline().strip().split()
    x_min, x_max = float(l[0]), float(l[1])
    
    l = file.readline().strip().split()
    y_min, y_max = float(l[0]), float(l[1])
    
    l = file.readline().strip().split()
    z_min, z_max = float(l[0]), float(l[1])

    if open:
        file.close()

    return (x_min, x_max, y_min, y_max, z_min, z_max)


def read_number_of_atoms(location: str) -> int:
    with open(location, "r") as f:
        i = 0
        for lines in f:
            i += 1
            if i == 4 :
                return int(lines.strip().split()[0])
            else:
                continue




def find_max(file: str, coord: str) -> float:
    
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


def find_min(file: str, coord: str) -> float:
    
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


def scale(target_volume: float, current_volume: float, num_walls: int) -> float:
    # calculate ratio between current and new box edges lengths
    if num_walls < 3:
        return pow(target_volume / current_volume, 1 / (3 - num_walls))
    if num_walls == 3:    
        return 1


class Atom:
    def __init__(self, id: int, x: float, y: float, z: float, type="A") -> None:
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
    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z


class Simulation_box:
    def __init__(self, x_min: float, x_max: float, y_min: float, y_max: float, z_min: float, z_max: float, atoms: int) -> None:
        self.min = Vector(x_min, y_min, z_min)
        self.max = Vector(x_max, y_max, z_max)
        self.x = x_max - x_min
        self.y = y_max - y_min
        self.z = z_max - z_min

        self.atoms = atoms
        self.mols = 1

        self.volume = self.x * self.y * self.z

class Molecule:
    def __init__(self, id: int, num_atoms: int) -> None:
        self.id = id   
        self.comp = []
        self.atoms = num_atoms

    def center(self) -> Atom:
        return self.comp[self.atoms//2]


    def center_of_mass(self) -> list:
        if len(self.comp) != self.atoms:
            raise Exception("Error: molecule not fully read." )

        com = np.zeros(3)

        for i in range(self.atoms):
            com += self.comp[i].position
            
        return com / self.atoms


    def polarization(self) -> None:
        return self.comp[self.atoms//2].position - (self.comp[-1].position + self.comp[0].position) / 2


    def shift(self, x: float, y: float, z: float) -> None:
        m = self.atoms // 2
        mid = Vector(self.comp[m].position[0], self.comp[m].position[1], self.comp[m].position[2] )
        for i in range(self.atoms):     
            self.comp[i].position[0] -= mid.x - x
            self.comp[i].position[1] -= mid.y - y
            self.comp[i].position[2] -= mid.z - z


    def rotate_x(self, theta: float) -> None:
        # change degrees to radians
        theta *= ( np.pi / 180 )
        for i in range(self.atoms):
            temp_y = self.comp[i].position[1]
            temp_z = self.comp[i].position[2]
            
            # self.comp[i].x = self.comp[i].x 
            self.comp[i].position[1] = temp_y * np.cos(theta) - temp_z * np.sin(theta)
            self.comp[i].position[2] = temp_y * np.sin(theta) + temp_z * np.cos(theta)


    def rotate_y(self, theta: float) -> None:
        # change degrees to radians
        theta *= ( np.pi / 180 )
        for i in range(self.atoms):
            temp_x = self.comp[i].position[0]
            temp_z = self.comp[i].position[2]

            self.comp[i].position[0] = temp_x * np.cos(theta) + temp_z * np.sin(theta)
            # self.comp[i].y = self.comp[i].y 
            self.comp[i].position[2] = -1 * temp_x * np.sin(theta) + temp_z * np.cos(theta)

    def rotate_z(self, phi: float) -> None:        
        # change degrees to radians
        phi *= ( np.pi / 180 )
        for i in range(self.atoms):
            temp_x = self.comp[i].position[0]
            temp_y = self.comp[i].position[1]
            
            self.comp[i].position[0] = temp_x * np.cos(phi) - temp_y * np.sin(phi)
            self.comp[i].position[1] = temp_x * np.sin(phi) + temp_y * np.cos(phi)
            # self.comp[i].z = self.comp[i].z 


    def print(self) -> None:
        for i in range(Molecule.atoms):
            print(f"{self.comp[i].position}")
        

    def director(self) -> list:
        director = self.comp[-1].position - self.comp[0].position
        director /= np.linalg.norm(director)
        
        return director


    def add(self, atom: Atom) -> None:
        self.comp.append(atom)

class Rescale:
    def __init__(self, current_packing_fraction: float, target_packing_fraction: float) -> None:
        self.current_packing_fraction = current_packing_fraction
        self.target_packing_fraction = target_packing_fraction

    target_box = Vector(1, 1, 1)
    target_volume = 1
    current_volume = 1
    factor = 1


class Line:
    def __init__(self, id: int, type: str, x: float, y: float, z: float) -> None:
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
    def __init__(self) -> None:
        self.components = []

    def assign(self, atom: Atom) -> None:
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
    def colour(self) -> float:
        return len(self.components)


class DirectorPixel(Pixel):
    def __init__(self) -> None:
        super().__init__()

    def Q(self) -> list[list]:
        Q = np.zeros((3,3))
        
        for director in self.components:
            Q += 0.5 * (3*np.tensordot(director, director, axes=0) - np.identity(3))
        
        if len(self.components) > 0:
            Q /= len(self.components)

        return Q

    def local_director(self) -> list:
        if len(self.components) == 0:
            return np.zeros(3)
        
        eigenvalue, eigenvector = scipy.sparse.linalg.eigsh(self.Q(), k=1, which="LM")
        eigenvector = np.reshape(eigenvector, -1)

        if eigenvector[-1] < 0:
            eigenvector *= -1

        return eigenvector


    def P2Value(self, direc=[1,1,1], trunc=True) -> float:
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
    
    def colour(self) -> str:
        # return self.P2Value()
        return self.local_director()[1]
        # return f"{self.local_director()[0]} {self.local_director()[1]} {self.local_director()[2]}"
    

class Screen:
    def __init__(self, x: int, y: int, pixel_class: Pixel) -> None:
        self.x = x
        self.y = y
        self.screen = [[pixel_class() for j in range(x)] for i in range(y)]

    def __repr__(self) -> str:
        z = 0.5
        data = str()
        for i in range(self.y):
            for j in range(self.x):
                pos = f"{z} {i/self.y} {j/self.x} " # z,y,x coords of pixel, 
                n = f"{len(self.screen[i][j].components)} "
                colour = f"{self.screen[i][j].colour()} "
                pix = pos + n + colour + '\n'

                data = data + pix

        return data


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

    def assign(self, atom: Atom, x: float, y: float) -> None:
        self.screen[y][x].assign(atom)

    def append_screenshot(self, screenshot: list) -> None:
        for i in range(self.y):
            for j in range(self.x):
                self.screen[i][j].merge_pixels(screenshot.screen[i][j])
    
    def colour(self) -> list[list]:
        colour = np.zeros((self.x, self.y))

        for i in range(self.y):
            for j in range(self.x):
                colour[i, j] = self.screen[i][j].colour()

        return colour

    def avg_director(self) -> list:
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
    def __init__(self, x: int, y: int, pixel_class: Pixel) -> None:
        super().__init__(x, y, pixel_class)


    def pixels_to_scroll(self, pixel_num: int, box: Simulation_box, drift: float) -> int:
        # determine numer of pixels that the image should be scrolled in order to get smectic interferences at the same place
        return int(pixel_num * drift / box.z)
    
    
    def scroll(self, pixels_to_scroll: int, side="both") -> None:
        # scrolling as slicing screenshot into two pieces and changing their order 
        if side == 'r':
            middle = self.x//2

            # split screenshot into left and right halfs
            L, R = [], []
            for i in range(self.x):
                L.append(self.screen[i][:middle])
                R.append(self.screen[i][middle:])
            
            # scroll right-hand-side of the screenshot
            R = R[pixels_to_scroll:] + R[:pixels_to_scroll]

            for i in range(self.x):
                self.screen[i] = L[i] + R[i]

        elif side == 'l':
            middle = self.x//2

            # split screenshot into left and right halfs
            L, R = [], []
            for i in range(self.x):
                L.append(self.screen[i][:middle])
                R.append(self.screen[i][middle:])
            
            # scroll left-hand-side of the screenshot
            L = L[pixels_to_scroll:] + L[:pixels_to_scroll]

            for i in range(self.x):
                self.screen[i] = L[i] + R[i]

        elif side == "both":
            self.screen = self.screen[pixels_to_scroll:] + self.screen[:pixels_to_scroll]
        else:
            raise Exception(f"Side can be 'l' or 'r' or 'both'!")


class HorizontalSlice(Screen):
    # initializing Slice as a single pixel of type matching current screen
    def __init__(self, screen: Screen, row: int) -> None:
        self.component = type(screen.screen[0][0])()
        self.row = row
        self.max_row = screen.x

    def __repr__(self) -> str:
        z = 0.5
        data = str()

        pos = f"{self.row / self.max_row} " # normalized index of slice, 
        d = self.slice_director()
        director = f"{d[0]} {d[1]} {d[2]}"
        data = data + pos + director + '\n'
        
        return data

    def read_slice(self, screen: Screen, start: int, end: int) -> None:
        for i in range(start, end):
            self.component.merge_pixels(screen.screen[self.row][i])
        
    def slice_director(self) -> list:
        if not isinstance(self.component, DirectorPixel):
            raise Exception("Pixel has to be Director Pixel")

        return self.component.local_director()


class VerticalSlice(Screen):
    # initializing Slice as a single pixel of type matching current screen
    def __init__(self, screen: Screen, row: int) -> None:
        self.component = type(screen.screen[0][0])()
        self.row = row

    def read_slice(self, screen: Screen, start: int, end: int) -> None:
        for i in range(start, end):
            self.component.merge_pixels(screen.screen[i][self.row])
        
    def slice_director(self) -> list:
        if not isinstance(self.component, DirectorPixel):
            raise Exception("Pixel has to be Director Pixel")

        return self.component.local_director()


class SmecticParameter():
    def __init__(self, periods: int) -> None:
        self.parameter = 0 + 0j
        self.count = 0
        self.periods = periods

    def __repr__(self) -> str:
       return f"{np.abs(self.parameter):.3f} | {self.count}"

    # only for vertical slices in y-z plane!
    def add_atom(self, center_coords: list, box: Simulation_box) -> None:
        self.parameter += np.exp(self.periods * 2*np.pi*1j * center_coords[-1] / box.z)
        self.count += 1
        # print(np.exp(2*np.pi*1j * center[-1] / box.z))

    def normalize(self) -> None:
        if self.count != 0:
            self.parameter /= self.count
    
    def read_screen(self, screen: Screen, start: int, end: int) -> None:
        HARD_CODED_LIMIT = 0
        # HARD_CODED_LIMIT = int((self.periods % 1.0) * screen.y / self.periods)
        # print(HARD_CODED_LIMIT, screen.y)

        for x in range(start, end):
        # for x in range(start, start+1):
            for y in range(HARD_CODED_LIMIT, screen.y):
                # for coords in screen.screen[y][x].components:
                #     self.add_atom(coords, box)
                self.parameter += np.exp(self.periods * 2*np.pi*1j * y / screen.y) * screen.screen[y][x].colour()
                self.count += screen.screen[y][x].colour()
                # print(self)

def read_matrix(file: str, n_slices: float, n_periods: int) -> list:
    list_of_params = np.zeros(n_slices, dtype=np.complex128)
    list_of_counts = np.zeros(n_slices)

    with open(file, "r") as f:
        for pixel in f:
            # read colour and coordinates of pixels
            z, y, x, N = pixel.split()
            z, y, x, N = float(z), float(y), float(x), int(N)

            # calculate slice to which pixel belongs and accumulate its contribution to smectic parameter
            index = int(x * n_slices)
            list_of_params[index] += np.exp(n_periods* 2*np.pi*1j * y) * N
            list_of_counts[index] += N

    # normalize parameters, take care of zero division
    for i in range(n_slices):
        if list_of_counts[i] != 0:
            list_of_params[i] /= list_of_counts[i]

    return list_of_params


def wrap_atom_to_box(atom: Atom, box: Simulation_box) -> Atom:
    atom.position[0] = atom.position[0] % box.x
    atom.position[1] = atom.position[1] % box.y
    atom.position[2] = atom.position[2] % box.z

    return atom

        

