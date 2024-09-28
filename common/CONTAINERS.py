import os
import numpy as np
from mmap import mmap

class Data:
    def __init__(self, location) -> None:
        self.location = location
        self.size = self.get_size()

    def get_location(self):
        return self.location
    
    def get_size(self):
        file_stats = os.stat(self.location)
        file_size = file_stats.st_size
        return file_size

class DataChunk(Data):
    def __init__(self, location, chunk_ID) -> None:
        super().__init__(location)
        self.ID = chunk_ID
        self.offset = self.ID * self.get_size()
    
    def get_ID(self):
        return self.ID

    def get_offset(self):
        return self.ID * self.get_size()
    
    def get_size(self):
        global NP
        chunk_size = self.get_size() // NP  
        chunk_size += mmap.ALLOCATIONGRANULARITY
        chunk_size -= chunk_size % mmap.ALLOCATIONGRANULARITY
        return chunk_size
    

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

class PhaseTracker:
    def __init__(self) -> None:
        C_left = 0 + 0j
        C_right = 0 + 0j
        C = 0 + 0j