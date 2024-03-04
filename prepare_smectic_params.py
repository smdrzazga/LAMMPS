import numpy as np
import matplotlib.pyplot as plt
import os
import banana_lib as sz
from scipy.optimize import curve_fit


desktop = r"C:\Users\Szymek\Desktop"
paths = {
    "desktop": desktop,
    "matrices": desktop + r"\LAMMPS_matrices",
    "matrices_centers": desktop + r"\LAMMPS_matrices\centers_matrices",
    "matrices_directors": desktop + r"\LAMMPS_matrices\directors_matrices",
    "smectic_params": desktop + r"\smectic_params_new.txt",
    "images": desktop + r"\checkpoint_01.03.24\images",
    "example": desktop + r"\LAMMPS_matrices\centers_matrices\centers_screen_bulk_6k_0.32.txt",
    "test": desktop + r"\test.png"
    }

N_SLICES = 30
DIRECTOR_PERIODS = 1
SMECTIC_PERIODS = 4
BOX_X = 70/8.62


def create_plot(location: str, target: str) -> None:
    # get packing fraction from file name 
    if type(location) == os.DirEntry: 
        location = location.name

    density = location.split("_")[-1].split(".")
    density = f"{density[0]}.{density[1]}"

    # create lists for plotting local smectic parameters
    smectic_params = sz.read_matrix(location, N_SLICES, SMECTIC_PERIODS)

    # data for plotting in box size units of length
    # print(smectic_params)
    y_axis = np.abs(smectic_params)
    x_axis = [i / N_SLICES for i in range(N_SLICES//2)]
    x_axis.extend([i / N_SLICES for i in reversed(range(N_SLICES//2))])

    # curve fitting
    # popt, pvar = curve_fit(lambda x, A, L: A*np.exp(-x / L), x_axis, y_axis, sigma=[0.005]*len(x_axis), absolute_sigma=True)
    # popt[1] *= box.x / x     # rescale pixels to atom diameters
    popt, pvar = curve_fit(lambda x, A, L: A*np.exp(-x / L), x_axis, y_axis, absolute_sigma=False)
    print(f"A = {popt[0]} \nlambda = {popt[1]}")
    print(f"var_A = {np.sqrt(pvar[0][0])} \nvar_lambda = {np.sqrt(pvar[1][1])}")
    
    with open(target, "a") as t:
        print(f"d = {density}", file=t)
        print(f"A = {popt[0]} \nvarA = {np.sqrt(pvar[0][0])} \nlambda = {popt[1]} \nvarLambda = {np.sqrt(pvar[1][1])} \n", file=t)

    # exponential plot
    plt.figure()
    plt.scatter(x_axis , y_axis, c='black')
    plt.xlabel(r"$\mathrm{distance\;to\;wall\;\,[box\;\,width]}$", fontsize=20)
    plt.ylabel(r"$\tau$", fontsize=24)
    plt.title(r'smectic order parameter, d\:=\:{}'.format(density))
    # plt.savefig(f"smectic_order_parameter_{density}.jpg")
    # plt.show()



if __name__ == "__main__":
    with open(paths["smectic_params"], "w+") as t: print('', end='', file=t)   # clear file

    for file in os.scandir(paths["matrices_centers"]):
        file = paths["matrices_centers"] + '\\' + file.name
        print(file)
        create_plot(file, paths["smectic_params"])
