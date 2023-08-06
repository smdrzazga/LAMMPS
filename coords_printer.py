import banana_lib as sz
import numpy as np

# arc-banana-shaped molecule angles and radius
molecule = sz.Molecule(1, 11)
chi = 130
alpha = (180 - chi)*np.pi/180
R = 1 / (2 * np.sin(alpha / 2 / molecule.atoms))

for i in range(molecule.atoms):
    x = R * np.cos( (i - molecule.atoms // 2) * alpha / molecule.atoms) - R * np.sin(alpha)
    y = 1
    z = R * np.sin( (i - molecule.atoms // 2) * alpha / molecule.atoms) + R * np.cos(alpha)

    atom = sz.Atom(i+1, 1, x, y, z)
    molecule.comp.append(atom)

# box needs to be inflated in order to achieve target packing fraction
# initial size of the simulation box
x = 140
y = 34
z = 45
packingFractionEnd = 0.308

# computing final box volume from target packing fraction and volume of all molecules
grid = sz.Vector(100, 30, 4)
volBoxStart = x * y * z
volAtoms =  grid.x*grid.y*grid.z * molecule.atoms * sz.Atom.volume
packingFractionStart = volAtoms / volBoxStart
volBoxEnd = volAtoms / packingFractionEnd 
scale = sz.scale( volBoxEnd, volBoxStart, 1 )
final_box = sz.Vector(x, y*scale, z*scale)

# variables determining lattice, on which molecules will be placed and distances between them
offset_mult = sz.Vector((final_box.x-8)/grid.x, (final_box.y-2)/grid.y, (z-2)/grid.z)
offset_add = sz.Vector(3, 1, -1)

# position of middle atom of molecule
mid = molecule.atoms // 2
position = sz.Vector( molecule.comp[mid].position[0], molecule.comp[mid].position[1], molecule.comp[mid].position[2] ) 

# translating molecule to be centered at (0,0,0), then rotating and translating back to original position
molecule.shift(0, 0, 0)
molecule.rotate_x(20)
molecule.shift(position.x, position.y, position.z)

# heading for LAMMPS read_file function
target = "C:/Users/Szymek/Desktop/tworzenie_atomow_z_pliku.txt"
sz.write_heading(target, molecule.atoms * grid.x*grid.y*grid.z, final_box)

# loops for printing atoms one at a time
id = 0
with open(target, "a") as f:
    for k in range(grid.z):
        # rotate whole layer of molecules by 90 deg for ever layer in z axis to get periodic structure every 45 vertical distance units
        position = sz.Vector( molecule.comp[mid].position[0], molecule.comp[mid].position[1], molecule.comp[mid].position[2] ) 
        molecule.shift(0,0,0)
        molecule.rotate_z(90)
        molecule.shift(position.x, position.y, position.z)
        
        for i in range(grid.x):
            for j in range(grid.y):
                for a in range(molecule.atoms):
                    id += 1
                    mol_id = i * grid.y * grid.z + j * grid.z + k + 1
                    x = molecule.comp[a].position[0] + offset_mult.x * i + offset_add.x
                    y = molecule.comp[a].position[1] + offset_mult.y * j + offset_add.y
                    z = molecule.comp[a].position[2] + offset_mult.z * k + offset_add.z
                    print("%d %d 1 %.5f %.2f %.5f" % (id, mol_id, x, y, z), file=f)

print(f"Packing fraction end: {packingFractionEnd} \nAll done!")
