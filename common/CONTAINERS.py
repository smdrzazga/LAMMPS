from MOLECULES import *
from mmap import mmap
import os

class SimulationData:
    def __init__(self, location) -> None:
        self.location = location
        self.size = self.get_file_size()

    def get_location(self):
        return self.location
    
    def get_file_size(self):
        file_stats = os.stat(self.location)
        file_size = file_stats.st_size
        return file_size

class File:
    def __init__(self, location) -> None:
        self.data = SimulationData(location)
        self.file = None

    def __del__(self):
        if self.file is not None:
            self.close()

    def open(self):
        open(self.data.get_location(), 'r+')

    def get_location(self):
        return self.data.location

    def get_file_size(self):
        return self.data.get_file_size()

    def close(self):
        self.file.close()
        self.file = None

    def clear(self):
        with open(self.file, 'r+w') as f:
            f.write('')

class ChunkData(SimulationData):
    def __init__(self, file_params: map, chunk_ID) -> None:
        super().__init__(file_params["location"])
        self.params = file_params
        self.ID = chunk_ID
        self.offset = self.get_offset()
    
    def get_ID(self):
        return self.ID

    def get_offset(self):
        start_percent, stop_percent = self.params["ANALYZE_RANGE"]
        return self.get_file_size()*start_percent + self.get_chunk_size() * self.ID
    
    def get_chunk_size(self):
        file_percent = self.params["ANALYZE_RANGE"][1] - self.params["ANALYZE_RANGE"][0]
        chunk_size = self.get_file_size() * file_percent // self.params["NP"]  
        return self._round_up(chunk_size, mmap.ALLOCATIONGRANULARITY)
    
    def _round_up(self, value, precision):
        value += precision
        value -= value % precision
        return value
    

class SimulationBox:
    def __init__(self, boundaries: list, n_atoms) -> None:
        self.update_boundaries(boundaries)

        self.volume = np.prod(self.size)
        self.atoms = int(n_atoms)

    def get_all_side_lengths(self) -> list:
        return self.size

    def get_side_length(self, i) -> float:
        return self.size[i]
    
    def get_num_atoms(self):
        return self.atoms
    
    def get_volume(self):
        return self.volume
    
    def update_boundaries(self, new_boundaries: list):
        [x_min, x_max, y_min, y_max, z_min, z_max] = new_boundaries
        self.min = np.array([x_min, y_min, z_min], dtype=np.float32)
        self.max = np.array([x_max, y_max, z_max], dtype=np.float32)
        self.size = np.array([x_max - x_min, y_max - y_min, z_max - z_min], dtype=np.float32)

