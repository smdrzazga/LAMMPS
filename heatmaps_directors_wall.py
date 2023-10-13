import numpy as np
import banana_lib as sz
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from multiprocessing.pool import Pool
import mmap
import time


# locations = [
#     "E:\LAMMPS_data/6k/all_snapshots_0.32.lammpstrj",
#     "E:\LAMMPS_data/6k/all_snapshots_0.318.lammpstrj",    
#     "E:\LAMMPS_data/6k/all_snapshots_0.316.lammpstrj",
#     "E:\LAMMPS_data/6k/all_snapshots_0.315.lammpstrj",
#     "E:\LAMMPS_data/6k/all_snapshots_0.314.lammpstrj",
#     "E:\LAMMPS_data/6k/all_snapshots_0.312.lammpstrj",
#     "E:\LAMMPS_data/6k/all_snapshots_0.31.lammpstrj",
#     "E:\LAMMPS_data/6k/all_snapshots_0.308.lammpstrj",
#     "E:\LAMMPS_data/6k/all_snapshots_0.305.lammpstrj",
#     "E:\LAMMPS_data/6k/all_snapshots_0.3.lammpstrj",
#     "E:\LAMMPS_data/6k/all_snapshots_0.29.lammpstrj",
#     "E:\LAMMPS_data/6k/all_snapshots_0.28.lammpstrj"
# ]

# locations = ["C:/Users/komok/Desktop/two_domain_tester/mirror_domains/two_domain_continuation/all_snapshots_0.32.lammpstrj"]

locations = [
    "E:/LAMMPS_data/double_z/all_snapshots_0.312.lammpstrj",
    "E:/LAMMPS_data/double_z/all_snapshots_0.32.lammpstrj"
]
NP = 6
SIZE = mmap.ALLOCATIONGRANULARITY * 2000
# input bin size and side of simulation box
x = 150
z = 150
plane = "xz"
DIRECTOR_PERIODS = 2

colorbar_limits = (-1, 1)


def analyze_batch(n, location):
    screen = sz.Screen(x, z, sz.DirectorPixel)
    screenshotDirector = sz.Screenshot(x, z, sz.DirectorPixel)
    screenshotCenter = sz.Screenshot(x, z, sz.CenterPixel)
    box = sz.Simulation_box(*sz.read_boundaries(location), sz.read_number_of_atoms(location))

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
                
                # if center.position[0] < 5:
                director = sz.Atom(molecule.id, *molecule.director())

                # assign director to the bin corresponding to the position of middle atom of the molecule
                pixel_position = screen.determine_pixel(center, box, plane)
                screenshotDirector.assign(director, *pixel_position)
                screenshotCenter.assign(center, *pixel_position)


            # if there is only one molecule remaining to read the full snapshot then execute following
            if atom.id == box.atoms:
                # eliminate Goldstone's mods by shifting whole system along z axis by:  L_z * Arg(C) / 2pi
                flow_left = box.z / DIRECTOR_PERIODS * np.angle(C_left) / (2*np.pi)
                flow_right = box.z / DIRECTOR_PERIODS * np.angle(C_right) / (2*np.pi)
                flow = box.z / DIRECTOR_PERIODS * np.angle(C) / (2*np.pi)
                pix_to_scroll_left = screenshotDirector.pixels_to_scroll(z, box, flow_left)
                pix_to_scroll_right = screenshotDirector.pixels_to_scroll(z, box, flow_right)
                pix_to_scroll_both = screenshotDirector.pixels_to_scroll(z, box, flow)
                
                screenshotDirector.scroll(pix_to_scroll_left, side="l")
                screenshotDirector.scroll(pix_to_scroll_right, side="r")
                # screenshotDirector.scroll(pix_to_scroll_both, side="both")

                # add corrected screenshot to the final image
                screen.append_screenshot(screenshotDirector)


                # check whether average director aligns with z direction
                # print(screenshotDirector.avg_director())
                # print(screen.avg_director())

                # clear current screenshot and flow measuring number C
                screenshotDirector = sz.Screenshot(x, z, sz.DirectorPixel)
                screenshotCenter = sz.Screenshot(x, z, sz.CenterPixel)
                C_left = 0 + 0j
                C_right = 0 + 0j
                C = 0 + 0j

                # print(pix_to_scroll_left, pix_to_scroll_right)
                # print(pix_to_scroll_both)

    print(f"Task is done!: {n}")
    return screen

def create_heatmap(screen: sz.Screen, location):
    # instantiate box for getting its size
    box = sz.Simulation_box(*sz.read_boundaries(location), sz.read_number_of_atoms(location))

    # create heatmap with a color scale
    plt.figure()
    if plane == "yz":
        plt.imshow(screen.colour(), interpolation="bilinear", clim = colorbar_limits, extent=(0, box.y, 0, box.z), aspect="auto")
    elif plane == "xz":
        plt.imshow(screen.colour(), interpolation="bilinear", clim = colorbar_limits, extent=(0, box.x, 0, box.z), aspect="auto")
    # if iter == 0:
    plt.colorbar()

    # plot labels and title
    plt.xlabel(f"{plane[0]}  [atom diameters]")
    plt.ylabel(f"{plane[1]}  [atom diameters]")

    density = location.split("_")[-1].split(".")
    zeros = ""
    for i in range(3-len(density[1])):
        zeros += "0"
    density = density[0] + "." + density[1] + zeros
    plt.title(f"directors, wall, d={density}")

    # show OR save plot - not both at once!
    # plt.show()
    plt.savefig(f"directors_{density}_total.png")



if __name__ == '__main__':

    for location in locations:
        screen = sz.Screen(x, z, sz.DirectorPixel)
        n = 90
        t1 = time.time()
        with Pool(NP) as executor:
            for result in executor.starmap(analyze_batch, zip([i for i in range(n)], [location]*n)):
            # for result in executor.starmap(analyze_batch, [(locations[0], i) for i in range(n)]):
                screen.append_screenshot(result)

        t2 = time.time()
        print(f"Time elapsed: {t2 - t1}")
        print("Here comes the heatmap!")
        create_heatmap(screen, location)
        print("Bye bye, heatmap!")
