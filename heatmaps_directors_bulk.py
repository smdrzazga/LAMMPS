from common.SCREEN import * 
from common.IO import *
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from multiprocessing.pool import Pool
import time


# locations = ["G:/lammps dane/npt/d_0.310/p_0.95_3x/all_snapshots_0.31.lammpstrj"]
# locations = ["C:/Users/Szymek/Desktop/middle_snapshot_5000000.lammpstrj"]
# locations = ["C:/Users/Szymek/Desktop/praca magisterska/kod/nematyk/all_snapshots_0.32.lammpstrj"]
locations = ["G:/lammps dane/two_domains/all_snapshots_0.32.lammpstrj"]


NP = 10
BATCH_START = 3
BATCH_STOP = 10
DIRECTOR_PERIODS = 2

AT_WALL = False
plane = "xz"
size = [150, 150]
ATOMS_IN_MOLECULE = 11

class DirectorScreenContainer:
    def __init__(self, size: tuple) -> None:
        self.screen = Screen(size, DirectorPixel)
        self.screenshotDirector = Screenshot(size, DirectorPixel)
        self.screenshotCenter = Screenshot(size, CenterPixel)
        self.molecule = Banana(1, ATOMS_IN_MOLECULE)


class BatchAnalyzer:
    def __init__(self, location) -> None:
        self.location = location
        self.screen_container = DirectorScreenContainer(size)

    def analyze_batch(self, ID):
        data_chunk = DataChunk(self.location, ID)
        reader = LAMMPSReader(data_chunk)
        phase_tracker = PhaseTracker()

        with open(location, "r+") as f:
            data = mmap.mmap(f.fileno(), length = SIZE, offset = (n)*SIZE)

            boundaries = sz.read_boundaries(data, N_ATOMS + 10, open=False, is_mmap=True)
            box = sz.SimulationBox(*boundaries, N_ATOMS)

            while True:
                line = data.readline().decode()

                if not line:
                    break

                # read atoms one by one from file 
                try:
                    atom = sz.Atom( line.split()[0], *line.split()[-3:])
                except:
                    continue

                # create molecule every 11 atoms
                if atom.id % 11 == 1:
                    molecule = sz.Molecule(atom.id//11 + 1, 11)
                
                molecule.add(atom)

                # if molecule is fully read then
                if len(molecule.comp) == molecule.atoms:
                    # translate molecule center back to simulation box 
                    center = sz.Atom(atom.id, *molecule.center_of_mass())
                    center = sz.wrap_atom_to_box(center, box)
            
                    # calculate C = sum_i p_y(i) exp (2 n pi z(i)/L_z)\
                    C += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)
                    if center.position[0] < box.x//2:
                        C_left += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)
                    else:
                        C_right += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)

                    # reject if center is not close to the wall, else add to screenshot
                    if not AT_WALL or (AT_WALL and center.position[0] < 5):
                        # assign director to the bin corresponding to the position of middle atom of the molecule
                        director = sz.Atom(molecule.id, *molecule.director())
                        pixel_position = screen.determine_pixel(center, box, plane)

                        screenshotDirector.assign(director, *pixel_position)
                        screenshotCenter.assign(center, *pixel_position)


                # if there is only one molecule remaining to read the full snapshot then execute following
                if atom.id == box.atoms:
                    # eliminate Goldstone's mods by shifting whole system along z axis by:  L_z * Arg(C) / 2pi
                    # flow_left = box.z / DIRECTOR_PERIODS * np.angle(C_left) / (2*np.pi)
                    # flow_right = box.z / DIRECTOR_PERIODS * np.angle(C_right) / (2*np.pi)
                    # pix_to_scroll_left = screenshotCenter.pixels_to_scroll(z, box, flow_left)
                    # pix_to_scroll_right = screenshotCenter.pixels_to_scroll(z, box, flow_right)
                    # screenshotDirector.scroll(pix_to_scroll_left, side="l")
                    # screenshotDirector.scroll(pix_to_scroll_right, side="r")

                    # flow = box.z / DIRECTOR_PERIODS * np.angle(C) / (2*np.pi)
                    # pix_to_scroll_both = screenshotCenter.pixels_to_scroll(z, box, flow)
                    # screenshotDirector.scroll(pix_to_scroll_both, side="both")

                    # add corrected screenshot to the final image
                    screen.append_screenshot(screenshotDirector)

                    # check whether average director aligns with z direction
                    # print(screenshotDirector.avg_director())
                    # print(screen.avg_director())

                    # clear current screenshot, box and flow measuring number C
                    boundaries = sz.read_boundaries(data, N_ATOMS + 10, open=False, is_mmap=True)
                    box = sz.SimulationBox(*boundaries, N_ATOMS)
                    screenshotDirector = sz.Screenshot(x, z, sz.DirectorPixel)
                    screenshotCenter = sz.Screenshot(x, z, sz.CenterPixel)
                    C, C_left, C_right = 0. + 0.j, 0. + 0.j, 0. + 0.j

        print(f"Task is done!: {n}")
        return screen
    




def create_scatter(screen: sz.Screen, location: str, start: int = 4, end: int = 9, new_file: bool = True) -> None:
    # common choices are 150x150 pix screen size
    # pixels 4:9 for wall and 73:78 for bulk  
    
    # initialize
    slices = [sz.HorizontalSlice(screen, row) for row in range(0, z)]
    slice_directors = np.zeros((z, 3))


    for i, slice in enumerate(slices):
        slice.read_slice(screen, start, end)
        slice_directors[i] = slice.component.local_director()

    # create empty file / clear existing values
    if new_file:
        with open(location, "w+") as l:
            print('', end='', file=l)

    with open(location, "a+") as l:
        # print(slice, end='', file=l)
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