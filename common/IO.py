from CONTAINERS import *
from SCREEN import Screen


class File(FileData):
    def __init__(self, location) -> None:
        super().__init__(location)
        self.file = None

    def __del__(self):
        if self.file is not None:
            self.close()

    def open(self):
        open(self.location(), 'r+')

    def close(self):
        self.file.close()
        self.file = None

    def clear(self):
        with open(self.file, 'r+w') as f:
            f.write('')


class Reader(File):
    def __init__(self, location) -> None:
        super().__init__(location)

    def open(self):
        f = open(self.location(), 'r+')
        self.file = mmap.mmap(f.fileno(), length=self.get_size())

    def get_line(self):
        try:
            return self.file.readline().decode()
        except:
            raise ValueError("File is not opened!")
        
    def get_line_split(self):
        return self.get_line().split()


class LAMMPSReader(Reader):
    def __init__(self, location) -> None:
        super().__init__(location)

    def open(self, batch_ID):
        chunk = DataChunk(self.location, batch_ID)
        f = open(self.get_location(), 'r+')
        self.file = mmap.mmap(f.fileno(), length=chunk.get_size(), offset=chunk.get_offset())

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

        return np.array([x_min, x_max, y_min, y_max, z_min, z_max], dtype=np.float32)

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


class Writer(File):
    def __init__(self, target_location) -> None:
        super().__init__(target_location)

    def open(self):
        self.file = open(self.location(), 'a+')
        
    def write_line(self, line):
        try:
            self.file.write(line)
        except:
            raise ValueError("File is not opened!")
    
    def print_obj(self, object):
        try:
            print(object, end='', file=self.file)
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

    