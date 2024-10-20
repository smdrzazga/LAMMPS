from common.CONTAINERS import *
from common.SCREEN import Screen
from config import ProcessingParameters
from mmap import mmap


class EmptyFile:
    def __init__(self, location) -> None:
        self.location = location
        self.file = None

    def open(self):
        self.file = open(self.get_location(), 'r+')
        return self.file 

    def get_location(self):
        return self.location

    def close(self):
        self.file.close()
        self.file = None

    def clear(self):
        with open(self.file, 'w') as f:
            f.write('')


class File(EmptyFile):
    def __init__(self, location) -> None:
        super().__init__(location)
        self.data = SimulationData(location)

    def get_file_size(self):
        return self.data.get_file_size()

    def close(self):
        self.file.close()
        self.file = None



class Reader:
    def __init__(self, location) -> None:
        self.file = File(location)
        self.map = None
        self.f = None

    def open(self):
        self.f = self.file.open()
        self.map = mmap(self.f.fileno(), length=self.file.get_file_size())

    def close(self):
        try:
            self.map.close()
            self.f.close()
        except:
            raise FileNotFoundError("File is already closed.")

    def get_line(self):
        try:
            return self.map.readline()
        except:
            raise ValueError("File is not opened!")
        
    def get_line_split(self):
        return self.get_line().split()
    
    def get_location(self):
        return self.file.get_location()


class LAMMPSReader(Reader):
    def __init__(self, proc_params: ProcessingParameters) -> None:
        super().__init__(proc_params.params["INPUT_FILE"])
        self.chunk = ChunkData(proc_params)

    def open(self, batch_ID):
        f = self.file.open()
        self.map = mmap(f.fileno(), length=self.chunk.get_chunk_size(), offset=self.chunk.get_offset(batch_ID))

    def read_boundaries(self) -> list:      
        l = self.find_line(self.is_box_header)
        x_min, x_max = np.array(l, dtype=np.float32)
        
        l = self.get_line_split()
        y_min, y_max = np.array(l, dtype=np.float32)
        
        l = self.get_line_split()
        z_min, z_max = np.array(l, dtype=np.float32)

        return np.array([x_min, x_max, y_min, y_max, z_min, z_max], dtype=np.float32)

    def read_number_of_atoms(self):
        atom_number_line = self.find_line(self.is_atoms_number_header)
        return int(atom_number_line[0])
    
    def read_atom_elements(self):
        atom_elements = self.find_line(self.is_atom_elements_header)
        return atom_elements[2:]

    def find_line(self, critera_func):
        line_split = self.get_line_split()
        while not critera_func(line_split):
            line_split = self.get_line_split()
        return self.get_line_split()

    def is_atom_elements_header(self, line_split):
        if len(line_split) < 5:
            return False
        return line_split[1] == 'ATOMS'

    def is_box_header(self, line_split):
        if len(line_split) < 2:
            return False
        return line_split[1] == 'BOX'
                
    def is_atoms_number_header(self, line_split):
        if len(line_split) < 4:
            return False
        return line_split[3] == 'ATOMS'


class Writer:
    def __init__(self, target_location) -> None:
        self.file = EmptyFile(target_location)
        self.pointer = None

    def open(self, mode='w+'):
        self.pointer = open(self.file.get_location(), mode)
        
    def close(self):
        self.pointer.close()

    def write(self, line, end=''):
        try:
            print(line, end=end, file=self.pointer)
        except:
            raise ValueError("File is not opened!")


class ScreenPrinter:
    def __init__(self, screen: Screen) -> None:
        self.screen = screen

    def print_screen(self, target_location):
        writer = Writer(target_location)
        writer.clear()
        writer.open()
        writer.print_obj(self.screen)
        writer.close()



# TODO: factorize
class LAMMPSWriter(Writer):
    def __init__(self) -> None:
        super().__init__()

    def write_heading(target: str, n: int, trim) -> None:
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

    