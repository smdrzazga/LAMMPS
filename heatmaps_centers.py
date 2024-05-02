import numpy as np
import banana_lib as sz
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from multiprocessing.pool import Pool
import mmap
import time

# locations = [
#     "G:/lammps dane/6k/all_snapshots_0.318.lammpstrj",    
#     "G:/lammps dane/6k/all_snapshots_0.316.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.315.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.314.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.312.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.31.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.3.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.305.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.29.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.28.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.308.lammpstrj",
#     "G:/lammps dane/6k/all_snapshots_0.32.lammpstrj"
# ]
locations = ["G:/lammps dane/4z_local/4z_240/all_snapshots_0.3.lammpstrj"]
# locations = [
#     "G:/lammps dane/double_z/all_snapshots_0.32.lammpstrj",
#     "G:/lammps dane/double_z/all_snapshots_0.312.lammpstrj"
# ]


NP = 10
SIZE = mmap.ALLOCATIONGRANULARITY * 2000
# input bin size and side of simulation box
x = 150
z = 150
plane = "xz"
DIRECTOR_PERIODS = 1
N_BATCH = 40

# colorbar_limits = (500, 1100)
# colorbar_limits = (1400, 1800)
# colorbar_limits = (30, 70)
# colorbar_limits = (80, 140)
colorbar_limits = (None, None)


def analyze_batch(n, location):
    box = sz.Simulation_box(*sz.read_boundaries(location), sz.read_number_of_atoms(location))
    screen = sz.Screen(x, z, sz.CenterPixel)
    screenshot = sz.Screenshot(x, z, sz.CenterPixel)

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

                # reject if center is not close to the wall, else add to screenshot
                # if molecule.center_of_mass()[0] < 5:
                screenshot.assign(center, *screenshot.determine_pixel(center, box, plane))    


            # if there is only one molecule remaining to read the full snapshot then execute following
            if atom.id == box.atoms:
                # eliminate Goldstone's mods by shifting whole system along z axis by:  L_z * Arg(C) / 2pi
                flow = box.z / DIRECTOR_PERIODS * np.angle(C) / (2*np.pi)
                pix_to_scroll_both = screenshot.pixels_to_scroll(z, box, flow)
                # screenshot.scroll(pix_to_scroll_both, side="both")
                
                # flow_left = box.z / DIRECTOR_PERIODS * np.angle(C_left) / (2*np.pi)
                # flow_right = box.z / DIRECTOR_PERIODS * np.angle(C_right) / (2*np.pi)
                # pix_to_scroll_left = screenshot.pixels_to_scroll(z, box, flow_left)
                # pix_to_scroll_right = screenshot.pixels_to_scroll(z, box, flow_right)
                # screenshot.scroll(pix_to_scroll_left, side="l")
                # screenshot.scroll(pix_to_scroll_right, side="r")
                

                # add corrected screenshot to the final image
                screen.append_screenshot(screenshot)

                # clear current screenshot and number C measuring flow 
                screenshot = sz.Screenshot(x, z, sz.CenterPixel)
                C_left = 0 + 0j
                C_right = 0 + 0j
                C = 0 + 0j



    print(f"Task is done!: {n}")
    return screen


def create_heatmap(screen: sz.Screen, location):
    # instantiate box for getting its size
    box = sz.Simulation_box(*sz.read_boundaries(location), sz.read_number_of_atoms(location))

    # plotting
    plt.figure()
    if plane == "yz":
        plt.imshow(screen.normalized_colour()[:, 4:-4], interpolation="bicubic", clim = colorbar_limits, extent=(0, box.y/sz.Molecule(1, 11).length(), 0, box.z/sz.Molecule(1, 11).length()), aspect="auto")
    elif plane == "xz":
        plt.imshow(screen.normalized_colour()[:, 4:-4], interpolation="bicubic", clim = colorbar_limits, extent=(0, box.x/sz.Molecule(1, 11).length(), 0, box.z/sz.Molecule(1, 11).length()), aspect="auto")        
    # if iter == 0:
    plt.colorbar()

    plt.xlabel(f"{plane[0]}  [molecule lengths]")
    plt.ylabel(f"{plane[1]}  [molecule lengths]")

    density = location.split("_")[-1].split(".")
    zeros = ""
    for i in range(3-len(density[1])):
        zeros += "0"
    density = density[0] + "." + density[1] + zeros
    step = location.split('_')[-1].split('.')[0]
    id = location.strip().split("/")[-2]
    # plt.title(f"centers, density={density}")

    # show OR save plot - not both at once!
    # plt.show()

    plt.savefig(f"centers_{density}.png")





if __name__ == '__main__':

    for location in locations:
        density = location.split('_')[-1].split('.')[0] + '.' + location.split('_')[-1].split('.')[1]
        mode = location.split('/')[-2]
        screen_file = "C:/Users/Szymek/Desktop/LAMMPS_matrices/centers_matrices/centers_screen_bulk_" + mode + '_' + density + ".txt"
        # screen_file = "C:/Users/Szymek/Desktop/test/centers_screen_bulk_" + mode + '_' + density + ".txt"

        screen = sz.Screen(x, z, sz.CenterPixel)
        t1 = time.time()
        with Pool(NP) as executor:
            i = 0
            for result in executor.starmap(analyze_batch, zip([i for i in range(N_BATCH)], [location]*N_BATCH)):
            # for result in executor.starmap(analyze_batch, [(locations[0], i) for i in range(N_BATCH)]):
                i += 1
                if i < 20: 
                    continue
                screen.append_screenshot(result)

        with open(screen_file, "w+") as t:
            print("", end='', file=t)

        with open(screen_file, "a+") as t:
            print(screen, end='', file=t)

        t2 = time.time()
        print(f"Time elapsed: {t2 - t1}")

        # print("Here comes the heatmap!")
        # create_heatmap(screen, location)
        # print("Bye bye, heatmap!")

