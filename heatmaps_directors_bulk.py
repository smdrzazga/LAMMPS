from common.SCREEN import * 
from common.IO import *
from common.UTILS import *
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
    def __init__(self, location, ID) -> None:
        self.location = location
        self.batch_data = DataChunk(location, ID)
        self.reader = LAMMPSReader(self.batch_data)
        self.binner = AtomBinner(size, view_plane=plane)
        self.phase_tracker = NTBPhaseTracker(DIRECTOR_PERIODS)

        self.screen = Screen(size, DirectorPixel)
        self.screenshotDirector = Screenshot(size, DirectorPixel)
        self.screenshotCenter = Screenshot(size, CenterPixel)
    
    def clear_screenshots(self):
        self.screenshotDirector = Screenshot(size, DirectorPixel)
        self.screenshotCenter = Screenshot(size, CenterPixel)

    def analyze_batch(self) -> Screen:
        self.reader.open()
        self.setup()

        line = self.reader.get_line_split()
        while line:
            try:
                atom = Atom(line[-3:], id=line[0])
            except:
                continue
            
            if self.molecule.is_full():
                self.molecule.clear()
            
            self.molecule.add(atom)

            if self.molecule.is_full:
                self.phase_tracker.update(self.molecule, self.box)
                center = Atom(self.molecule.center_of_mass())
                center.wrap(self.box.get_all_side_lengths())

                if ANALYSE_WALL and self.molecule.is_at_wall():
                    director = Atom(self.molecule.director())
                    pixel_position = self.binner.determine_pixel(center)

                    self.screenshotDirector.assign(director, pixel_position)
                    self.screenshotCenter.assign(center, pixel_position)

# ----------------------------------- CONTINUE HERE -------------------------------------------------------------

            if atom.is_last(N_ATOMS):
                if ELIMINATE_GOLDSTONE:
                    self._eliminate_goldstone_mods_both_sides()
                self.finalize_screenshot_analysis()
                
            line = self.reader.get_line_split()

        print(f"Task is done!: {self.batch_data.ID}")
        return screen

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

    def setup(self):
        boundaries = self.reader.read_boundaries()
        N_ATOMS = self.reader.read_number_of_atoms()
        self.box = SimulationBox(boundaries, N_ATOMS)
        self.molecule = Banana(1, ATOMS_IN_MOLECULE)

    def update_box(self):
        boundaries = self.reader.read_boundaries()
        self.box.update_boundaries(boundaries)

    def finalize_screenshot_analysis(self):
        self.screen.append_screenshot(self.screenshotDirector)
        self.update_box()
        self.clear_screenshots()
        self.phase_tracker.clear()


def create_scatter(screen: Screen, location: str, start: int = 4, end: int = 9, new_file: bool = True) -> None:
    slices = [sz.HorizontalSlice(screen, row) for row in range(0, z)]
    slice_directors = np.zeros((z, 3))


    for i, slice in enumerate(slices):
        slice.read_slice(screen, start, end)
        slice_directors[i] = slice.component.local_director()

    if new_file:
        with open(location, "w+") as l:
            print('', end='', file=l)

    with open(location, "a+") as l:
        for i in range(z):
            print(f"{i} {slice_directors[i, 0]} {slice_directors[i, 1]} {slice_directors[i, 2]}", file=l)



if __name__ == '__main__':

    for location in locations:
        N_ATOMS = sz.read_number_of_atoms(location)
        N_BATCH = BATCH_STOP - BATCH_START
        screen = sz.Screen(x, z, sz.DirectorPixel)

        density = location.split('_')[-1].split('.')[0] + '.' + location.split('_')[-1].split('.')[1]
        mode = location.split('/')[-2]
        screen_file = "C:/Users/Szymek/Desktop/LAMMPS_matrices/directors_screen_bulk_" + mode + '_' + density + ".txt"

        t1 = time.time()
        with Pool(NP) as executor:
            args = zip([i for i in range(BATCH_START, BATCH_STOP)], [location]*N_BATCH, [N_ATOMS]*N_BATCH)
            for result in executor.starmap(analyze_batch, args):
                screen.append_screenshot(result)

        # with open(screen_file, "w+") as t:
        #     print("", end='', file=t)

        # with open(screen_file, "a+") as t:
        #     print("Here comes the heatmap!")
        #     print(screen, end='', file=t)
        #     print("Finished! Yay!")

        # t2 = time.time()
        # print(f"Time elapsed: {t2 - t1}")
        # print("Bye bye, heatmap!")

    # screen_file = "C:/Users/Szymek/Desktop/3d_0.305.txt"
    # create_scatter(screen, screen_file, start=4, end=9)
    # for i in range(10, 145, 5):
    #     start = i
    #     end = i + 5
    #     print("Here comes the heatmap!")

    #     create_scatter(screen, screen_file, start=start, end=end, new_file=False)
    #     print("Bye bye, heatmap!")

        for s in range (20, 40):
            window = 20
            e = s + window

            screen_file = "C:/Users/Szymek/Desktop/scatters/" + str(s) + "-" + str(e) + "_0.32.txt"
            create_scatter(screen, screen_file, start=s, end=e)