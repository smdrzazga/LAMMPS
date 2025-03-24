import numpy as np
import banana_lib as sz
from multiprocessing.pool import Pool
import mmap
import time
import os


# locations = ["G:/lammps dane/new_chi/p1.05/all_snapshots_1.05.lammpstrj"]
locations = ["G:/lammps dane/two_domains/all_snapshots_2x.lammpstrj"]
# locations = ["C:/Users/Szymek/Desktop/ref_snap.lammpstrj"]



NP = 11
# input data and side of simulation box
BATCH_START = 20
BATCH_STOP = 90
DIRECTOR_PERIODS = 1
SIZE = mmap.ALLOCATIONGRANULARITY * 2000

AT_WALL = False
plane = "xz"
x = 150
z = 150


def analyze_batch(n, location, N_ATOMS):
    screen = sz.Screen(x, z, sz.DirectorPixel)
    screenshotDirector = sz.Screenshot(x, z, sz.DirectorPixel)
    screenshotCenter = sz.Screenshot(x, z, sz.CenterPixel)
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
                atom = sz.Atom( line.split()[0], [line.split()[-3:]])
            except:
                continue

            # create molecule every 11 atoms
            if atom.id % 11 == 1:
                molecule = sz.Molecule(atom.id//11 + 1, 11)
            
            molecule.add(atom)

            # if molecule is fully read then
            if len(molecule.comp) == molecule.atoms:
                # translate molecule center back to simulation box 
                center = sz.Atom(atom.id, molecule.center_of_mass())
                center = sz.wrap_atom_to_box(center, box)
        
                # calculate C = sum_i p_y(i) exp (2 n pi z(i)/L_z)\
                C += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)
                if center.position[0] < box.x//2:
                    C_left += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)
                else:
                    C_right += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)

                # reject if center is not close to the wall, else add to screenshot
                if not AT_WALL or (AT_WALL and center.position[0] > 115):
                    # assign director to the bin corresponding to the position of middle atom of the molecule
                    director = sz.Atom(molecule.id, molecule.director())
                    pixel_position = screen.determine_pixel(center, box, plane)
    
                    screenshotDirector.assign(director, *pixel_position)
                    screenshotCenter.assign(center, *pixel_position)


            # if there is only one molecule remaining to read the full snapshot then execute following
            if atom.id == 3*box.atoms//2:
                # eliminate Goldstone's mods by shifting whole system along z axis by:  L_z * Arg(C) / 2pi
                # flow_left = box.z / DIRECTOR_PERIODS * np.angle(C_left) / (2*np.pi)
                # flow_right = box.z / DIRECTOR_PERIODS * np.angle(C_right) / (2*np.pi)
                # pix_to_scroll_left = screenshotCenter.pixels_to_scroll(z, box, flow_left)
                # pix_to_scroll_right = screenshotCenter.pixels_to_scroll(z, box, flow_right)
                # screenshotDirector.scroll(pix_to_scroll_left, side="l")
                # screenshotDirector.scroll(pix_to_scroll_right, side="r")

                flow = box.z / DIRECTOR_PERIODS * np.angle(C) / (2*np.pi)
                pix_to_scroll_both = screenshotCenter.pixels_to_scroll(z, box, flow)
                # screenshotDirector.scroll(pix_to_scroll_both, side="both")

                # add corrected screenshot to the final image
                screen.append_screenshot(screenshotDirector)

                # check whether average director aligns with z direction
                # print(screenshotDirector.avg_director())
                # print(screen.avg_director())

                # clear current screenshot, box and flow measuring number C
                boundaries = sz.read_boundaries(data, N_ATOMS + 10, open=False, is_mmap=True)
                box = sz.Simulation_box(*boundaries, N_ATOMS)
                screenshotDirector = sz.Screenshot(x, z, sz.DirectorPixel)
                screenshotCenter = sz.Screenshot(x, z, sz.CenterPixel)
                C, C_left, C_right = 0. + 0.j, 0. + 0.j, 0. + 0.j

    print(f"Task is done!: {n}")
    return screen


def create_director_matrix(screen: sz.Screen, location):
    results = []
    with Pool(NP) as executor:
        args = zip([(i, j) for i in range(z) for j in range(x)],  [screen.screen[i][j] for i in range(z) for j in range(x)])
        for result in executor.starmap(_get_director, args):
            results.append(result)

    with open(location, "w+") as l:
        for line in results:
            i, j, director = line
            print(f"{i} {j} {director[0]} {director[1]} {director[2]}", file=l)
                
def _get_director(coords, pixel):
    i, j = coords
    return [i, j, pixel.local_director()]



if __name__ == '__main__':

    for location in locations:
        # N_ATOMS = sz.read_number_of_atoms(location)
        N_ATOMS = 165000
        N_BATCH = BATCH_STOP - BATCH_START
        screen = sz.Screen(x, z, sz.DirectorPixel)

        density = location.split('_')[-1].split('.')[0] + '.' + location.split('_')[-1].split('.')[1]
        mode = location.split('/')[-2]
        screen_file = "C:/Users/" + os.getlogin() + "/Desktop/LAMMPS_matrices/directors_matrices/directors_" + ("second_wall" if AT_WALL else "bulk") + "_" + mode + '_' + density + ".txt"

        t1 = time.time()
        with Pool(NP) as executor:
            args = zip([i for i in range(BATCH_START, BATCH_STOP)], [location]*N_BATCH, [N_ATOMS]*N_BATCH)
            for result in executor.starmap(analyze_batch, args):
                screen.append_screenshot(result)

        print("Here comes the heatmap!")
        create_director_matrix(screen, screen_file)
        print("Finished! Yay!")

        t2 = time.time()
        print(f"Time elapsed: {t2 - t1}")
        print("Bye bye, heatmap!")
