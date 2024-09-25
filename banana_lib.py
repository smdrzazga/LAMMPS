from math import sqrt
import os
import numpy as np
import functools
from scipy.sparse.linalg import eigsh
import scipy.sparse.linalg
from mmap import mmap



def scale(target_volume: float, current_volume: float, num_walls: int) -> float:
    # calculate ratio between current and new box edges lengths
    if num_walls < 3:
        return pow(target_volume / current_volume, 1 / (3 - num_walls))
    if num_walls == 3:    
        return 1

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



class DataChunk:
    chunk_size = 0
    def __init__(self, location, nth_Chunk) -> None:
        self.location = location
        self.ID = nth_Chunk

    def get_offset(self):
        return self.ID * self.get_size()
    
    def get_ID(self):
        return self.ID
    
    def get_size(self):
        global NP        
        file_stats = os.stat(self.location)
        file_size = file_stats.st_size
        return file_size // NP
    
    def get_location(self):
        return self.location
    
    

class BatchReader:
    def __init__(self, data_chunk: DataChunk) -> None:
        self.chunk = data_chunk
        self.file = None

    def __del__(self):
        if self.file is not None:
            self.close()

    def open(self):
        f = open(self.chunk.get_location(), 'r+')
        self.file = mmap.mmap(f.fileno(), length=self.chunk.get_size(), offset=self.chunk.get_offset())

    def close(self):
        self.file.close()

    def get_line(self):
        return self.file.readline().decode()

    def get_line_split(self):
        return self.get_line().split()

    def read_boundaries(self) -> list:
        line = self.get_line()
                        
        while not self.is_box_header(line):
            line = self.get_line()
        
        l = self.get_line_split()
        x_min, x_max = np.array(l, dtype=np.float32)
        
        l = self.get_line_split()
        y_min, y_max = np.array(l, dtype=np.float32)
        
        l = self.get_line_split()
        z_min, z_max = np.array(l, dtype=np.float32)

        return [x_min, x_max, y_min, y_max, z_min, z_max]

    def read_number_of_atoms(self):
        self.open()
        line_split = self.get_line_split()
        
        while not self.is_atoms_number_header(line_split):
            line_split = self.get_line_split()
        atom_number_line = self.get_line_split()

        return int(atom_number_line.strip()[0])

    def is_box_header(self, line_split):
        if len(line_split) < 2:
            return False
        return line_split[1] == 'BOX'
                
    def is_atoms_number_header(self, line_split):
        if len(line_split) < 4:
            return False
        return line_split[3] == 'ATOMS'


class SimulationBox:
    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max, n_atoms) -> None:
        self.min = np.array([x_min, y_min, z_min], dtype=np.float32)
        self.max = np.array([x_max, y_max, z_max], dtype=np.float32)
        self.size = np.array([x_max - x_min, y_max - y_min, z_max - z_min], dtype=np.float32)

        self.atoms = int(n_atoms)
        self.mols = 1

        self.volume = np.prod(self.size)

    def get_all_side_lengths(self) -> list:
        return self.size

    def get_side_length(self, i) -> float:
        return self.size[i]


class Atom:
    volume = 1.333 * 3.14 * 0.5 * 0.5 * 0.5

    def __init__(self, id: int, position: list, type=1) -> None:
        self.id = int(id)
        self.type = type
        self.position = np.array(position, dtype=np.float32)

    def __repr__(self) -> str:
        return f"{self.id} | {self.type} | {self.position}"


class Molecule:
    def __init__(self, id: int, num_atoms: int) -> None:
        self.id = id   
        self.comp = []
        self.atoms = num_atoms

    def add(self, atom: Atom) -> None:
        self.comp.append(atom)

    def center_atom(self) -> Atom:
        return self.comp[self.atoms//2]

    def center_of_mass(self) -> list:
        com = np.zeros(3)
        for i in range(self.atoms):
            com += self.comp[i].position
        return com / self.atoms

    def is_molecule_full(self) -> bool:
        return len(self.comp) == self.atoms

    def set_position(self, position: list) -> None:
        mid = self.center_of_mass()
        for i in range(len(self.comp)):     
            self.comp[i].position += position - mid

    def shift(self, displacement: list) -> None:
        for i in range(len(self.comp)):     
            self.comp[i].position += displacement

    def rotate_x(self, theta_deg: float) -> None:
        theta = theta_deg * np.pi/180
        for i in range(self.atoms):
            temp_y = self.comp[i].position[1]
            temp_z = self.comp[i].position[2]
            
            self.comp[i].position[1] = temp_y * np.cos(theta) - temp_z * np.sin(theta)
            self.comp[i].position[2] = temp_y * np.sin(theta) + temp_z * np.cos(theta)

    def rotate_y(self, theta_deg: float) -> None:
        theta = theta_deg * np.pi/180
        for i in range(self.atoms):
            temp_x = self.comp[i].position[0]
            temp_z = self.comp[i].position[2]

            self.comp[i].position[0] = temp_x * np.cos(theta) + temp_z * np.sin(theta)
            self.comp[i].position[2] = -1 * temp_x * np.sin(theta) + temp_z * np.cos(theta)

    def rotate_z(self, phi: float) -> None:        
        phi *= ( np.pi / 180 )
        for i in range(self.atoms):
            temp_x = self.comp[i].position[0]
            temp_y = self.comp[i].position[1]
            
            self.comp[i].position[0] = temp_x * np.cos(phi) - temp_y * np.sin(phi)
            self.comp[i].position[1] = temp_x * np.sin(phi) + temp_y * np.cos(phi)     


class Banana(Molecule):
    def __init__(self, id: int, num_atoms: int) -> None:
        super().__init__(id, num_atoms)

    def director(self) -> list:
        director = self.comp[-1].position - self.comp[0].position
        director /= np.linalg.norm(director)
        return director

    def polarization(self) -> None:
        middle_atom = self.comp[self.atoms//2].position
        average = (self.comp[-1].position + self.comp[0].position) / 2
        return middle_atom - average


class Rescale:
    def __init__(self, current_packing_fraction: float, target_packing_fraction: float) -> None:
        self.current_packing_fraction = current_packing_fraction
        self.target_packing_fraction = target_packing_fraction

    target_box = np.array([1, 1, 1])
    target_volume = 1
    current_volume = 1
    factor = 1


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
    def __init__(self, size: list, pixel_class: Pixel) -> None:
        self._assert_size_dimension(size, 2)
        self.size = np.array(size)
        self.screen = [[pixel_class() for j in range(size[0])] for i in range(size[1])]

    def _assert_size_dimension(self, size, expected):
        if len(size) != expected:
            return Exception("Screen size should have 2 values: horizontal & vertical!")

    def __repr__(self) -> str:
        z = 0.5
        data = str()
        for i in range(self.y):
            for j in range(self.x):
                pos = f"{z} {i/self.size[1]} {j/self.size[0]} " # z,y,x coords of pixel, 
                n = f"{len(self.screen[i][j].components)} "
                colour = f"{self.screen[i][j].colour()} "
                pix = pos + n + colour + '\n'

                data = data + pix

        return data

    def get_size_in_pixels(self):
        return self.size        

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


    def pixels_to_scroll(self, pixel_num: int, box: SimulationBox, drift: float) -> int:
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


class AtomBinner:
    _indice_map = {'x': 0, 'y': 1, 'z': 2}
    plane = ''

    def __init__(self, sim_box: SimulationBox, screen_size_in_pixels: list, view_plane='xz') -> None:
        self.plane = view_plane
        self.box = sim_box
        self.screen_size_in_pixels = np.array([screen_size_in_pixels])

    def determine_pixel(self, atom: Atom):
        box_size = self._get_box_size_from_viewing_plane()
        atom_position = self._get_atom_position_from_viewing_plane(atom)

        pixel_size = self._compute_pixel_size(box_size)
        pixel_coords = self._translate_atom_to_pixel(atom_position, pixel_size)
        pixel_coords = self._assign_outer_atoms_to_borders(pixel_coords)

        return pixel_coords
    

    def _get_box_size_from_viewing_plane(self):
        i, j = self._plane_to_coord_indices(self.plane)
        box_size = np.array([self.box.get_side_length(i), self.box.get_side_length(j)])

        return box_size

    def _get_atom_position_from_viewing_plane(self, atom: Atom):
        i, j = self._plane_to_coord_indices(self.plane)
        atom_position = np.array([atom.position[i], atom.position[j]])

        return atom_position
    
    def _plane_to_coord_indices(self, plane: str):        
        i = self._indice_map[plane[0]]
        j = self._indice_map[plane[1]]

        return i, j

    def _compute_pixel_size(self, box_size):
        pixel_size = box_size / self.screen_size_in_pixels

        return pixel_size

    def _translate_atom_to_pixel(self, atom_position, pixel_size):
        pixel_coords = np.array(atom_position // pixel_size, dtype=np.int32)

        return pixel_coords

    def _assign_outer_atoms_to_borders(self, pixel_coords):
        print(f"DEBUG: pixel coords before: {pixel_coords}\n")
        for i, pixel_coord in enumerate(pixel_coords):
            if pixel_coord >= self.screen_size_in_pixels[i]:
                pixel_coords[i] = self.screen_size_in_pixels[i] - 1
            elif pixel_coord < 0:
                pixel_coords[i] = 0
            else:
                pass

        print(f"DEBUG: pixel coords after: {pixel_coords}\n")
        return pixel_coords


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
    
    def get_local_director(self):
        return self.component.local_director()


class SmecticParameter():
    def __init__(self, periods: int) -> None:
        self.parameter = 0 + 0j
        self.count = 0
        self.periods = periods

    def __repr__(self) -> str:
       return f"{np.abs(self.parameter):.3f} | {self.count}"

    # only for vertical slices in y-z plane!
    def add_atom(self, center_coords: list, box: SimulationBox) -> None:
        self.parameter += np.exp(self.periods * 2*np.pi*1j * center_coords[-1] / box.z)
        self.count += 1

    def normalize(self) -> None:
        if self.count != 0:
            self.parameter /= self.count
    
    def read_screen(self, screen: Screen, start: int, end: int) -> None:
        HARD_CODED_LIMIT = 0

        for x in range(start, end):
            for y in range(HARD_CODED_LIMIT, screen.y):
                self.parameter += np.exp(self.periods * 2*np.pi*1j * y / screen.y) * screen.screen[y][x].colour()
                self.count += screen.screen[y][x].colour()


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


def wrap_atom_to_box(atom: Atom, box: SimulationBox) -> Atom:
    atom.position[0] = atom.position[0] % box.x
    atom.position[1] = atom.position[1] % box.y
    atom.position[2] = atom.position[2] % box.z

    return atom

        

