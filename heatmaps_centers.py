import numpy as np
import banana_lib as sz
import matplotlib.pyplot as plt

location = "C:/Users/komok/Desktop/all_snapshots_0.32.lammpstrj"
# location = "C:/Users/komok/Desktop/all_snapshots_0.31.lammpstrj"
# location = "C:/Users/komok/Desktop/all_snapshots_0.3.lammpstrj"
# location = "C:/Users/komok/Desktop/all_snapshots_0.29.lammpstrj"
# location = "C:/Users/komok/Desktop/programy/middle_snapshot_100000.lammpstrj"
target = "C:/Users/komok/Desktop/check.txt"

# input bin size and side of simulation box
x = 150
# y = 11
z = 150
plane = "yz"

screen = sz.Screen(x, z, sz.CenterPixel)
screenshot = sz.Screenshot(x, z, sz.CenterPixel)
box = sz.Simulation_box(*sz.read_boundaries(location), sz.read_number_of_atoms(location))

i = 0
C = 0 + 0j
with open(location, "r") as f:
    with open(target, "w") as t:
        for lines in f:
            # read atoms one by one from file 
            try:
                atom = sz.Atom( *lines.strip().split() )
                i += 1
            except:
                continue

            # create molecule every 11 atoms
            if atom.id % 11 == 1:
                molecule = sz.Molecule(atom.id//11 + 1, 11)
            
            molecule.add(atom)

            # if molecule is fully read then
            if atom.id % 11 == 0:
                # calculate C = sum_i p_y(i) exp (2 n pi z(i)/L_z)
                C += molecule.polarization()[1] * np.exp(2j*1*np.pi * molecule.center_of_mass()[2] / box.z)
                
                # reject if center is not close to the wall
                if molecule.center_of_mass()[0] < 5:
                    # translate molecule center back to simulation box 
                    center = sz.wrap_atom_to_box(molecule.center(), box)
                    screenshot.assign(center, *screenshot.determine_pixel(center, box, plane))    
        

            # if there is only one molecule remaining to read the full snapshot then execute following
            if atom.id == box.atoms:
                # eliminate Goldstone's mods by shifting whole system along z axis by:  L_z * Arg(C) / 2pi
                flow = box.z * np.angle(C) / (2*np.pi)
                pix_to_scroll = screenshot.pixels_to_scroll(z, box, flow)
                screenshot.scroll(pix_to_scroll)

                # add corrected screenshot to the final image
                screen.append_screenshot(screenshot)

                # clear current screenshot and number C measuring flow 
                screenshot = sz.Screenshot(x, z, sz.CenterPixel)
                C = 0

                print(pix_to_scroll)


            # if i > 5e5:
            #     break


# fig, (ax1, ax2) = plt.subplots(1, 2)

# heatmap1 = ax1.imshow(screen.colour())
# ax1.set_title('total')
# fig.colorbar(heatmap1, ax=ax1)


# heatmap2 = ax2.imshow(screenshot.colour())
# ax2.set_title('screenshot')
# fig.colorbar(heatmap2, ax=ax2)

# plt.show()



# print(screenshot)

plt.imshow(screen.colour())
plt.colorbar()

plt.xlabel(f"{plane[0]}")
plt.ylabel(f"{plane[1]}")
plt.title("centers heatmap")

plt.show()

