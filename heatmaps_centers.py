import numpy as np
import banana_lib as sz
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from multiprocessing.pool import Pool
import mmap
import time


locations = ["G:/lammps dane/two_domains/all_snapshots_1.01.lammpstrj"]

# locations = [
#     "G:/lammps dane/6k/all_snapshots_0.32.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.318.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.316.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.315.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.314.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.312.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.31.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.308.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.3.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.29.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.28.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.305.lammpstrj"
# ]


NP = 11
# input data and side of simulation box
BATCH_START = 20
BATCH_STOP = 30
DIRECTOR_PERIODS = 1
SIZE = mmap.ALLOCATIONGRANULARITY * 1000

AT_WALL = False
plane = "xz"
x = 150
z = 150


def analyze_batch(n, location, N_ATOMS):
    screenshot = sz.Screenshot(x, z, sz.CenterPixel)
    screen = sz.Screen(x, z, sz.CenterPixel)
    molecule = sz.Molecule(1, 11)

    C_left = 0 + 0j
    C_right = 0 + 0j
    C = 0 + 0j

    with open(location, "r+") as f:
        data = mmap.mmap(f.fileno(), length = SIZE, offset = (n)*SIZE)

        boundaries = sz.read_boundaries(data, N_ATOMS + 10, open=False, is_mmap=True)
        box = sz.Simulation_box(*boundaries, N_ATOMS)

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
                molecule.translate(box.min.x, box.min.y, box.min.z)
                center = sz.Atom(atom.id, *molecule.center_of_mass())
                center = sz.wrap_atom_to_box(center, box)
                
                # calculate C = sum_i p_y(i) exp (2 n pi z(i)/L_z)\
                C += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)
                if center.position[0] < box.x//2:
                    C_left += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)
                else:
                    C_right += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)

                # reject if center is not close to the wall, else add to screenshot
                if not AT_WALL or (AT_WALL and center.position[0] > 65):
                    screenshot.assign(center, *screenshot.determine_pixel(center, box, plane))    


            # if there is only one molecule remaining to read the full snapshot then execute following
            if atom.id == box.atoms:
                # eliminate Goldstone's mods by shifting whole system along z axis by:  L_z * Arg(C) / 2pi
                flow = box.z / DIRECTOR_PERIODS * np.angle(C) / (2*np.pi)
                pix_to_scroll_both = screenshot.pixels_to_scroll(z, box, flow)
                # screenshot.scroll(pix_to_scroll_both, side="both")

                # add corrected screenshot to the final image
                screen.append_screenshot(screenshot)

                # clear current screenshot, box and number C measuring flow 
                boundaries = sz.read_boundaries(data, N_ATOMS + 10, open=False, is_mmap=True)
                box = sz.Simulation_box(*boundaries, N_ATOMS)
                screenshot = sz.Screenshot(x, z, sz.CenterPixel)
                C_left = 0 + 0j
                C_right = 0 + 0j
                C = 0 + 0j

    print(f"Task is done!: {n}")
    return screen


def main():
    for location in locations:
        N_ATOMS = sz.read_number_of_atoms(location)
        N_BATCH = BATCH_STOP - BATCH_START
        screen = sz.Screen(x, z, sz.CenterPixel)

        density = location.split('_')[-1].split('.')[0] + '.' + location.split('_')[-1].split('.')[1]
        mode = location.split('/')[-2]
        screen_file = "C:/Users/Szymek/Desktop/LAMMPS_matrices/centers_matrices/centers_screen_" + ("wall" if AT_WALL else "bulk") + "_" + mode + '_' + density + ".txt"
        # screen_file = "C:/Users/Szymek/Desktop/LAMMPS_matrices/centers_matrices/centers_screen_bulk_two_domains_1.01.txt"

        t1 = time.time()
        with Pool(NP) as executor:
            args = zip([i for i in range(BATCH_START, BATCH_STOP)], [location]*N_BATCH, [N_ATOMS]*N_BATCH)
            for result in executor.starmap(analyze_batch, args):
                screen.append_screenshot(result)

        print("Here comes the heatmap!")

        with open(screen_file, "w+") as t:
            print("", end='', file=t)

        with open(screen_file, "a+") as t:
            print(screen, end='', file=t)

        t2 = time.time()
        print(f"Time elapsed: {t2 - t1}")
        print("Bye bye, heatmap!")


if __name__ == '__main__':
    main()