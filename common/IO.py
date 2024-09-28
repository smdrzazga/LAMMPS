from CONTAINERS import *


class File(Data):
    def __init__(self) -> None:
        self.file = None

    def __del__(self):
        if self.file is not None:
            self.close()

    def open(self):
        f = open(self.get_location(), 'r+')
        self.file = mmap.mmap(f.fileno(), length=self.get_size())

    def close(self):
        self.file.close()

    def clear(self):
        with open(self.file, 'r+w') as f:
            f.write('')
    

class LAMMPSReader(File):
    def __init__(self, data_chunk: DataChunk) -> None:
        super().__init__()
        self.chunk = data_chunk

    def open(self):
        f = open(self.get_location(), 'r+')
        self.file = mmap.mmap(f.fileno(), length=self.get_size(), offset=self.chunk.get_offset())

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


# TODO: factorize
class LAMMPSWriter(File):
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

    