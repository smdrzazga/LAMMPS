from CONTAINERS import *
from MOLECULES import *
from scipy.sparse.linalg import eigsh



class Pixel:
    def __init__(self) -> None:
        self.components = []

    def _get_components(self):
        return self.components
    
    def num_elements(self):
        return len(self.components)

    def assign(self, atom: Atom) -> None:
        self.components.append(atom.position)

    def mean(self) -> float:
        mean = 0
        if len(self.components) > 0:
            mean = np.mean(self.components, axis=0)
            mean /= len(self.components)
        return mean
    
    def colour(self) -> float:
        return 1

    def merge_pixels(self, new_pixel) -> None:
        self.components.extend(new_pixel.components)

    def print(self):
        for atom in self.components:
            print(atom)


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
        
        eigenvalue, eigenvector = eigsh(self.Q(), k=1, which="LM")
        eigenvector = np.reshape(eigenvector, -1)

        if eigenvector[-1] < 0:
            eigenvector *= -1

        return eigenvector


    def P2Value(self, direc=[1,1,1]) -> float:
        direction = np.array(direc) / np.linalg.norm(np.array(direc))

        if len(self.components) == 0:
            return 0
        
        p2value = 0.5 * (3 * ( np.dot(self.local_director(), direction) )**2 - 1)
        p2value = (p2value + 0.5) / 1.5 

        return p2value
    
    def colour(self) -> str:
        # return self.P2Value()
        return self.local_director()[1]
    

class Screen:
    def __init__(self, size: list, pixel_class: Pixel) -> None:
        self._assert_size_dimension(size, 2)
        self.size = np.array(size)
        self.screen = [[pixel_class() for j in range(size[0])] for i in range(size[1])]
        self.pixel_type = pixel_class

    def _assert_size_dimension(self, size, expected):
        if len(size) != expected:
            return Exception("Screen size should have 2 values: horizontal & vertical!")

    def get_size_in_pixels(self):
        return self.size        
    
    def _get_pixel(self, x, y):
        return self.screen[y][x]
    
    def get_row(self, i):
        return self.screen[i]
    
    def _get_pixel_type(self):
        return self.pixel_type

    def __repr__(self) -> str:
        z = 0.5
        data = str()
        x, y = self.get_size_in_pixels()
        for i in range(x):
            for j in range(y):
                pos = f"{z} {j/y} {i/x} "
                n = f"{self._get_pixel(i,j).num_elements()} "
                colour = f"{self._get_pixel(i,j).colour()} "
                pix = pos + n + colour + '\n'

                data += pix
        return data
    
    def get_slice(self, range: list[list]):
        slice = np.array(range)
        return self.screen[slice[1,0] : slice[1,1]][slice[0,0]:slice[0,1]]

    def assign(self, atom: Atom, x: float, y: float) -> None:
        self._get_pixel(x,y).assign(atom)

    def append_screenshot(self, screenshot: list) -> None:
        x, y = self.get_size_in_pixels()
        for i in range(x):
            for j in range(y):
                self._get_pixel(i,j).merge_pixels(screenshot._get_pixel(i,j))
    
    def colour(self) -> list[list]:
        colour = np.zeros(self.size)
        x, y = self.get_size_in_pixels()
        for i in range(x):
            for j in range(y):
                colour[i, j] = self._get_pixel(i,j).colour()

        return colour

    def avg_director(self) -> list:
        try: 
            self._compute_avg_director()
        except:
            raise TypeError("Pixel has to be an instance of Director Pixel")
    
    def _compute_avg_director(self) -> list:
        Q = np.zeros((3,3))
        x, y = self.get_size_in_pixels()
        for i in range(x):
            for j in range(y):
                Q += self._get_pixel(i,j).Q()

        Q /= (np.prod(self.size))**2
        _, avg_director = eigsh(Q, k=1, which="LA")
        avg_director = np.reshape(avg_director, -1)

        return avg_director


class Screenshot(Screen):
    def __init__(self, size: list, pixel_class: Pixel) -> None:
        super().__init__(size, pixel_class)

    def pixels_to_scroll(self, pixel_num: int, box: SimulationBox, drift: float) -> int:
        return (pixel_num * drift) // box.get_side_length(2)

    def scroll_right_side(self, pixels_to_scroll: int) -> None:
        L = self._get_slice_left()
        R = self._get_slice_right()
        R = self._scroll_list(R, pixels_to_scroll)

        x, _ = self.get_size_in_pixels()
        for i in range(x):
            self.screen[i] = L[i] + R[i]

    def scroll_left_side(self, pixels_to_scroll: int) -> None:
        middle = self.get_size_in_pixels()[0] // 2
        
        L = self.get_slice([[None, middle], [None, None]])
        R = self.get_slice([[middle, None], [None, None]])
        L = self._scroll_list(L, pixels_to_scroll)

        x, _ = self.get_size_in_pixels()
        for i in range(x):
            self.screen[i] = L[i] + R[i]
        
    def scroll_both_sides(self, pixels_to_scroll: int) -> None:
        lower_part = self.get_slice([[None, None], [pixels_to_scroll, None]])
        upper_part = self.get_slice([[None, None], [None, pixels_to_scroll]])
        self.screen = lower_part + upper_part
    
    def _get_slice_left(self):
        middle = self.get_size_in_pixels()[0] // 2
        return self.get_slice([[None, middle], [None, None]])

    def _get_slice_right(self):
        middle = self.get_size_in_pixels()[0] // 2
        return self.get_slice([[middle, None], [None, None]])
    
    def _scroll_list(self, list, pixels_to_scroll):
        return list[pixels_to_scroll:] + list[:pixels_to_scroll]
    


class AtomBinner:
    _indice_map = {'x': 0, 'y': 1, 'z': 2}
    plane = ''

    def __init__(self, sim_box: SimulationBox, screen_size_in_pixels: list, view_plane='xz') -> None:
        self.plane = view_plane
        self.box = sim_box
        self.screen_size_in_pixels = np.array(screen_size_in_pixels, dtype=np.float32)

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
        for i, pixel_coord in enumerate(pixel_coords):
            if pixel_coord >= self.screen_size_in_pixels[i]:
                pixel_coords[i] = self.screen_size_in_pixels[i] - 1
            elif pixel_coord < 0:
                pixel_coords[i] = 0
            else:
                pass

        return pixel_coords




