import banana_lib as sz
import numpy as np

# arc-banana-shaped molecule angles and radius
molecule = sz.Molecule(1, 11)
chi = 130
alpha = (180 - chi)*np.pi/180
R = 1 / (2 * np.sin(alpha / 2 / molecule.atoms))

for i in range(molecule.atoms):
    x = R * np.cos( (i - molecule.atoms // 2) * alpha / molecule.atoms) - R * np.sin(alpha)
    y = R * np.sin( (i - molecule.atoms // 2) * alpha / molecule.atoms) + R * np.cos(alpha)
    z = 1

    atom = sz.Atom(i+1, x, y, z, type=1)
    molecule.comp.append(atom)

# box needs to be inflated in order to achieve target packing fraction
# initial size of the simulation box
x = 70
y = 45
z = 105
packingFractionEnd = 0.31
N_WALLS = 1     # REMEMBER ABOUT SETTING SCALE OF FINAL_BOX
N_MOLS_PER_PERIOD = 4
# IF THERE IS A PROBLEM WITH COLLAPSED MOLECULES INTO A PLANE - CHECK ATOM CONSTRUCTOR

# computing final box volume from target packing fraction and volume of all molecules
grid = sz.Vector(60, 40, 8)
volBoxStart = x * y * z
volAtoms =  grid.x*grid.y*grid.z * molecule.atoms * sz.Atom.volume
packingFractionStart = volAtoms / volBoxStart
volBoxEnd = volAtoms / packingFractionEnd 
scale = sz.scale( volBoxEnd, volBoxStart, N_WALLS )
final_box = sz.Vector(x, y*scale, z*scale)

# variables determining lattice, on which molecules will be placed and distances between them
offset_mult = sz.Vector((final_box.x-6)/grid.x, (final_box.y-0)/grid.y, (final_box.z-0)/grid.z)
offset_add = sz.Vector(0.5, 1, 0)

# position of middle atom of molecule
mid = molecule.atoms // 2
position = sz.Vector( *molecule.comp[mid].position ) 

# NTB
# translating molecule to be centered at (0,0,0), then rotating and translating back to original position
molecule.shift(0, 0, 0)
molecule.rotate_x(70)
molecule.shift(position.x, position.y, position.z)

# heading for LAMMPS read_file function
target = "C:/Users/Szymek/Desktop/tworzenie_atomow_z_pliku.txt"
sz.write_heading(target, molecule.atoms * grid.x*grid.y*grid.z, final_box)

# loops for printing atoms one at a time
id = 0
with open(target, "a") as f:
    for k in range(grid.z):
        # rotate whole layer of molecules by 90 deg for ever layer in z axis to get periodic structure every 45 vertical distance units
        position = sz.Vector( *molecule.comp[mid].position ) 
        molecule.shift(0,0,0)
        molecule.rotate_z(360//N_MOLS_PER_PERIOD)
        molecule.shift(position.x, position.y, position.z)
        
        for j in range(grid.x):
            for i in range(grid.y):
                for a in range(molecule.atoms):
                    id += 1
                    mol_id = k * grid.x * grid.y + j * grid.y + i + 1
                    x = molecule.comp[a].position[0] + offset_mult.x * j + offset_add.x
                    y = molecule.comp[a].position[1] + offset_mult.y * i + offset_add.y
                    z = molecule.comp[a].position[2] + offset_mult.z * k + offset_add.z
                    print("%d %d 1 %.5f %.2f %.5f" % (id, mol_id, x, y, z), file=f)

print(f"Packing fraction end: {packingFractionEnd} \nAll done!")


# #NSB
# # translating molecule to be centered at (0,0,0), then rotating and translating back to original position
# SB_ANGLE = 30
# molecule.shift(0, 0, 0)
# molecule.rotate_x(90)
# molecule.rotate_z(90)
# molecule.rotate_x(-1*SB_ANGLE)

# molecule.shift(position.x, position.y, position.z)

# # heading for LAMMPS read_file function
# target = "C:/Users/Szymek/Desktop/tworzenie_atomow_z_pliku.txt"
# sz.write_heading(target, molecule.atoms * grid.x*grid.y*grid.z, final_box)

# # loops for printing atoms one at a time
# id = 0
# with open(target, "a") as f:
#     for k in range(grid.z):
#         position = sz.Vector( *molecule.comp[mid].position ) 
#         molecule.shift(0,0,0)
#         if k%6 == 0 or k%6 == 1:
#             molecule.rotate_x(SB_ANGLE)
#         elif k%6 == 2:
#             molecule.rotate_z(180)
#             molecule.rotate_x(2*SB_ANGLE)
#         elif k%6 == 3 or k%6 == 4:
#             molecule.rotate_x(-1*SB_ANGLE)
#         else:
#             molecule.rotate_z(180)
#             molecule.rotate_x(-2*SB_ANGLE)
        
#         if k%3 == 2:
#             position.z -= 1.2
#         elif k%3 != 2:
#             position.z += 0.5

#         # molecule.shift(position.x, position.y, position.z)
#         molecule.shift(position.x, position.y, position.z)
        
#         for j in range(grid.x):
#             for i in range(grid.y):
#                 for a in range(molecule.atoms):
#                     id += 1
#                     mol_id = k * grid.x * grid.y + j * grid.y + i + 1
#                     x = molecule.comp[a].position[0] + offset_mult.x * j + offset_add.x
#                     y = molecule.comp[a].position[1] + offset_mult.y * i + offset_add.y
#                     z = molecule.comp[a].position[2] + offset_mult.z * k + offset_add.z
#                     print("%d %d 1 %.5f %.2f %.5f" % (id, mol_id, x, y, z), file=f)

# print(f"Packing fraction end: {packingFractionEnd} \nAll done!")