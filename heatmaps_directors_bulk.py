import numpy as np
import banana_lib as sz
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from multiprocessing.pool import Pool
import mmap
import time



# locations = [
#     "G:/lammps dane/6k/all_snapshots_0.316.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.315.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.31.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.305.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.29.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.28.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.32.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.314.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.312.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.307.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.318.lammpstrj",    
#     "G:/lammps dane/6k/all_snapshots_0.3.lammpstrj"
# ]
locations = ["G:/lammps dane/4z_local/4z_240/all_snapshots_0.3.lammpstrj"]
# locations = [
#     "G:/lammps dane/double_z/all_snapshots_0.312.lammpstrj",
#     "G:/lammps dane/double_z/all_snapshots_0.32.lammpstrj",
#       "G:/lammps dane/double_z/test_const_z_90/all_snapshots_0.298.lammpstrj"
# ]


NP = 10
DIRECTOR_PERIODS = 1
N_BATCH = 40
# input bin size and side of simulation box
x = 150
z = 150
plane = "xz"
SIZE = mmap.ALLOCATIONGRANULARITY * 2000 * DIRECTOR_PERIODS


def analyze_batch(n, location):
    box = sz.Simulation_box(*sz.read_boundaries(location), sz.read_number_of_atoms(location))
    screen = sz.Screen(x, z, sz.DirectorPixel)
    screenshotDirector = sz.Screenshot(x, z, sz.DirectorPixel)
    screenshotCenter = sz.Screenshot(x, z, sz.CenterPixel)

    molecule = sz.Molecule(1, 11)
    C_left = 0 + 0j
    C_right = 0 + 0j
    C = 0 + 0j
    with open(location, "r+") as f:
        data = mmap.mmap(f.fileno(), length = SIZE, offset = (n)*SIZE)
        molecule = sz.Molecule(1, 11)
        while True:
            # atom = sz.Atom( *read_atom(line).T)
            line = data.readline().decode()

            if not line:
                break

            # read atoms one by one from file 
            try:
                # atom = sz.Atom( *line.strip().split() )
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
                # center = sz.Atom(atom.id, 'A', *molecule.center_of_mass())
                center = sz.Atom(atom.id, *molecule.center_of_mass())
                center = sz.wrap_atom_to_box(center, box)
        
                # calculate C = sum_i p_y(i) exp (2 n pi z(i)/L_z)\
                C += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)
                if center.position[0] < box.x//2:
                    C_left += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)
                else:
                    C_right += molecule.polarization()[1] * np.exp(2j*DIRECTOR_PERIODS*np.pi * center.position[2] / box.z)

                # if middle.position[0] < 5:
                director = sz.Atom(molecule.id, *molecule.director())

                # assign director to the bin corresponding to the position of middle atom of the molecule
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

                # clear current screenshot and flow measuring number C
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
        density = location.split('_')[-1].split('.')[0] + '.' + location.split('_')[-1].split('.')[1]
        mode = location.split('/')[-2]
        screen_file = "C:/Users/Szymek/Desktop/LAMMPS_matrices/directors_screen_bulk_" + mode + '_' + density + ".txt"

        screen = sz.Screen(x, z, sz.DirectorPixel)
        t1 = time.time()
        with Pool(NP) as executor:
            i = 0
            for result in executor.starmap(analyze_batch, zip([i for i in range(N_BATCH)], [location]*N_BATCH)):
            # for result in executor.starmap(analyze_batch, [(locations[0], i) for i in range(N_BATCH)]):
                i += 1
                if i < 20: 
                    continue
                screen.append_screenshot(result)

        # with open(screen_file, "w+") as t:
        #     print("", end='', file=t)

        # with open(screen_file, "a+") as t:
        #     print("Here comes the heatmap!")
        #     print(screen, end='', file=t)
        #     print("Finished! Yay!")


    screen_file = "C:/Users/Szymek/Desktop/3d_0.300.txt"
    create_scatter(screen, screen_file, start=4, end=9)
    for i in range(10, 145, 5):
        start = i
        end = i + 5
        print("Here comes the heatmap!")

        create_scatter(screen, screen_file, start=start, end=end, new_file=False)
        print("Bye bye, heatmap!")

        

# 39, 46