from common.SCREEN import * 
from common.IO import *
from common.PROCESSING import *
from common.PROCESSING import *
from concurrent.futures import ProcessPoolExecutor
from multiprocessing.pool import Pool
import matplotlib.pyplot as plt
import time


# locations = ["G:/lammps dane/npt/d_0.310/p_0.95_3x/all_snapshots_0.31.lammpstrj"]
# locations = ["C:/Users/Szymek/Desktop/middle_snapshot_5000000.lammpstrj"]
# locations = ["C:/Users/Szymek/Desktop/praca magisterska/kod/nematyk/all_snapshots_0.32.lammpstrj"]
locations = ["G:/lammps dane/two_domains/all_snapshots_0.32.lammpstrj"]


NP = 10
BATCH_START = 3
BATCH_STOP = 10
DIRECTOR_PERIODS = 2

ATOMS_IN_MOLECULE = 11
ANALYSE_WALL = False
ELIMINATE_GOLDSTONE = False
plane = "xz"
size = [150, 150]


class BatchAnalyzer:
    def __init__(self, location) -> None:
        self.location = location
        self.reader = LAMMPSReader(location)
        self.binner = AtomBinner(size, view_plane=plane)
        self.phase_tracker = NTBPhaseTracker(DIRECTOR_PERIODS)

        self.screen = Screen(size, DirectorPixel)
        self.screenshotDirector = Screenshot(size, DirectorPixel)
        self.screenshotCenter = Screenshot(size, CenterPixel)

    def analyze_batch(self, ID) -> Screen:
        self.reader.open()
        self.setup()

        line = self.reader.get_line_split()
        while line:
            try:
                self.read_atom(line)
            except:
                continue
                     
            self.molecule.add(self.atom)

            if not self.molecule.is_full():
                continue

            if ANALYSE_WALL and not self.molecule.is_at_wall():
                continue

            self.add_molecule_to_pixel()
            self.molecule.clear()

            if self.atom.is_last(self.box.get_num_atoms()):
                if ELIMINATE_GOLDSTONE:
                    self._eliminate_goldstone_mods_both_sides()
                self.finalize_screenshot_analysis()
                
            line = self.reader.get_line_split()

        print(f"Task is done!: {ID}")
        return screen

    def setup(self):
        boundaries = self.reader.read_boundaries()
        N_ATOMS = self.reader.read_number_of_atoms()
        self.box = SimulationBox(boundaries, N_ATOMS)
        self.molecule = Banana(1, ATOMS_IN_MOLECULE)
        self.atom = Atom([0, 0, 0])

    def read_atom(self, split_line):
        self.atom = Atom(split_line[-3:], id=split_line[0])

    def add_molecule_to_pixel(self):
        self.phase_tracker.update(self.molecule, self.box)

        director = Atom(self.molecule.director())
        center = Atom(self.molecule.center_of_mass())
        center.wrap(self.box.get_all_side_lengths())
        pixel_position = self.binner.determine_pixel(center)

        self.screenshotDirector.assign(director, pixel_position)
        self.screenshotCenter.assign(center, pixel_position)

    def _eliminate_goldstone_mods_both_sides(self, tracker: NTBPhaseTracker, box: SimulationBox) -> None:
        drift = tracker.compute_com_drift(box)
        pix_to_scroll_both = self.screenshotCenter.pixels_to_scroll(size[1], box, drift)
        self.screenshotDirector.scroll_both_sides(pix_to_scroll_both)

    def _eliminate_goldstone_mods_right_side(self, tracker: NTBPhaseTracker, box: SimulationBox) -> None:
        drift = tracker.compute_com_drift(box)
        pix_to_scroll_both = self.screenshotCenter.pixels_to_scroll(size[1], box, drift)
        self.screenshotDirector.scroll_right_side(pix_to_scroll_both)
    
    def _eliminate_goldstone_mods_left_side(self, tracker: NTBPhaseTracker, box: SimulationBox) -> None:
        drift = tracker.compute_com_drift(box)
        pix_to_scroll_both = self.screenshotCenter.pixels_to_scroll(size[1], box, drift)
        self.screenshotDirector.scroll_left_side(pix_to_scroll_both)

    def finalize_screenshot_analysis(self):
        self.screen.append_screenshot(self.screenshotDirector)
        self.update_box()
        self.clear_screenshots()
        self.phase_tracker.clear()
    
    def update_box(self):
        boundaries = self.reader.read_boundaries()
        self.box.update_boundaries(boundaries)
    
    def clear_screenshots(self):
        self.screenshotDirector = Screenshot(size, DirectorPixel)
        self.screenshotCenter = Screenshot(size, CenterPixel)


if __name__ == '__main__':
    for location in locations:
        screen = Screen(size, DirectorPixel)

        density = location.split('_')[-1].split('.')[0] + '.' + location.split('_')[-1].split('.')[1]
        mode = location.split('/')[-2]
        target_file = "C:/Users/Szymek/Desktop/LAMMPS_matrices/directors_screen_bulk_" + mode + '_' + density + ".txt"

        t1 = time.time()
        with Pool(NP) as executor:
            IDs = [ID for ID in range(NP)]
            analyze = BatchAnalyzer(location).analyze_batch
            for result in executor.starmap(analyze, IDs):
                screen.append_screenshot(result)


            print("Here comes the heatmap!")
            ScreenPrinter(screen).print_screen(target_file)
            print("Finished! Yay!")

        t2 = time.time()
        print(f"Time elapsed: {t2 - t1}")
        print("Bye bye, heatmap!")
