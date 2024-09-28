import numpy as np

class Atom:
    volume = 1.333 * 3.14 * 0.5 * 0.5 * 0.5

    def __init__(self, position: list, id=0, type=1) -> None:
        self.id = int(id)
        self.type = type
        self.position = np.array(position, dtype=np.float32)

    def __repr__(self) -> str:
        return f"{self.id} | {self.type} | {self.position}"
    
    def wrap(self, size: list):
        self.position %= np.array(size)


class Molecule:
    def __init__(self, id: int, max_atoms: int) -> None:
        self.id = id   
        self.comp = []
        self.max_atoms = max_atoms
        self.current_atoms = 0

    def _get_atom(self, i) -> Atom:
        return self.comp[i]

    def add(self, atom: Atom) -> None:
        self.comp.append(atom)
        self.current_atoms += 1

    def clear(self):
        self.comp = []
        self.current_atoms = 0

    def center_atom(self) -> Atom:
        return self._get_atom(self.current_atoms//2)

    def center_of_mass(self) -> list:
        com = np.zeros(3)
        for i in range(self.current_atoms):
            com += self._get_atom(i).position
        return com / self.current_atoms

    def is_full(self) -> bool:
        return self.current_atoms == self.max_atoms
    
    def is_at_wall(self) -> bool:
        return self.center_of_mass[0] < 5

    def set_position(self, position: list) -> None:
        mid = self.center_of_mass()
        for i in range(self.current_atoms):     
            self._get_atom(i).position += position - mid

    def shift(self, displacement: list) -> None:
        for i in range(self.current_atoms):     
            self._get_atom(i).position += displacement

    def rotate_x(self, theta_deg: float) -> None:
        theta = theta_deg * np.pi/180
        for i in range(self.current_atoms):
            temp_y = self._get_atom(i).position[1]
            temp_z = self._get_atom(i).position[2]
            
            self._get_atom(i).position[1] = temp_y * np.cos(theta) - temp_z * np.sin(theta)
            self._get_atom(i).position[2] = temp_y * np.sin(theta) + temp_z * np.cos(theta)

    def rotate_y(self, theta_deg: float) -> None:
        theta = theta_deg * np.pi/180
        for i in range(self.current_atoms):
            temp_x = self._get_atom(i).position[0]
            temp_z = self._get_atom(i).position[2]

            self._get_atom(i).position[0] = temp_x * np.cos(theta) + temp_z * np.sin(theta)
            self._get_atom(i).position[2] = -1 * temp_x * np.sin(theta) + temp_z * np.cos(theta)

    def rotate_z(self, phi: float) -> None:        
        phi *= np.pi / 180
        for i in range(self.current_atoms):
            temp_x = self._get_atom(i).position[0]
            temp_y = self._get_atom(i).position[1]
            
            self._get_atom(i).position[0] = temp_x * np.cos(phi) - temp_y * np.sin(phi)
            self._get_atom(i).position[1] = temp_x * np.sin(phi) + temp_y * np.cos(phi)     


class Banana(Molecule):
    def __init__(self, id: int, num_atoms: int) -> None:
        super().__init__(id, num_atoms)

    def director(self) -> list:
        director = self._get_atom(-1).position - self._get_atom(0).position
        director /= np.linalg.norm(director)
        return director

    def polarization(self) -> None:
        middle_atom = self._get_atom(self.current_atoms//2).position
        average = (self._get_atom(-1).position + self._get_atom(0).position) / 2
        return middle_atom - average