import numpy as np
import banana_lib as sz
import matplotlib.pyplot as plt


# location = "C:/Users/komok/Desktop/all_snapshots_0.32.lammpstrj"
# location = "C:/Users/komok/Desktop/all_snapshots_0.31.lammpstrj"
location = "C:/Users/komok/Desktop/all_snapshots_0.3.lammpstrj"
# location = "C:/Users/komok/Desktop/all_snapshots_0.29.lammpstrj"
# location = "C:/Users/komok/Desktop/programy/middle_snapshot_100000.lammpstrj"
target = "C:/Users/komok/Desktop/check.txt"


# input bin size and side of simulation box
x = 150
# y = 11
z = 150
plane = "xz"

screen = sz.Screen(x, z, sz.DirectorPixel)
screenshotDirector = sz.Screenshot(x, z, sz.DirectorPixel)
screenshotCenter = sz.Screenshot(x, z, sz.CenterPixel)
box = sz.Simulation_box(*sz.read_boundaries(location), sz.read_number_of_atoms(location))

molecule_len = 3
molecule_mid = 1
i = 0
C = 0 + 0j

with open(location, "r") as f:
    with open(target, "w") as t:
        for lines in f:
            i += 1
            # read atoms one by one from file 
            try:
                atom = sz.Atom( *lines.strip().split() )
            except:
                continue

            # if atom is not first, middle or last atom in molecule then do not read it at all 
            if atom.id % 11 != 1 and atom.id % 11 != 0 and atom.id % 11 != 6:
                continue

            # create molecule every 11 atoms
            if atom.id % 11 == 1:
                molecule = sz.Molecule(atom.id//11 + 1, molecule_len)
            
            molecule.add(atom)

            # if molecule is fully read then
            if atom.id % 11 == 0:
                # calculate C = sum_i p_y(i) exp (2 n pi z(i)/L_z)
                C += molecule.polarization()[1] * np.exp(2j*1*np.pi * molecule.center_of_mass()[2] / box.z)
                
                # if molecule.comp[molecule_mid].position[0] < 5:

                director = sz.Atom(molecule.id, "A", *molecule.director())

                # translate atom back to simulation box and reject if it is not center of a molecule
                middle = molecule.comp[molecule_mid]
                middle = sz.wrap_atom_to_box(middle, box)

                # assign director to the bin corresponding to the position of middle atom of the molecule
                # x, y = screen.determine_pixel(middle, box, plane)
                screenshotDirector.assign(director, *screen.determine_pixel(middle, box, plane))
                screenshotCenter.assign(middle, *screen.determine_pixel(middle, box, plane))


            # if there is only one molecule remaining to read the full snapshot then execute following
            if atom.id == box.atoms:
                # eliminate Goldstone's mods by shifting whole system along z axis by:  L_z * Arg(C) / 2pi
                flow = box.z * np.angle(C) / (2*np.pi)
                pix_to_scroll = screenshotCenter.pixels_to_scroll(z, box, flow)
                screenshotDirector.scroll(pix_to_scroll)

                # add corrected screenshot to the final image
                screen.append_screenshot(screenshotDirector)


                # check whether average director aligns with z direction
                # print(screenshotDirector.avg_director())
                # print(screen.avg_director())

                # clear current screenshot and flow measuring number C
                screenshotDirector = sz.Screenshot(x, z, sz.DirectorPixel)
                screenshotCenter = sz.Screenshot(x, z, sz.CenterPixel)
                C = 0

                print(pix_to_scroll)



            # if i > 2e7:
            #     break



plt.imshow(screen.colour())
plt.colorbar()

plt.xlabel(f"{plane[0]}")
plt.ylabel(f"{plane[1]}")
plt.title("directors heatmap")

plt.show()


