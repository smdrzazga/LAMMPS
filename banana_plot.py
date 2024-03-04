import numpy as np
import banana_lib as sz
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



plt.rcParams.update({'font.size': 18})
plt.rcParams.update({"text.usetex": True})


class Plotter():
    def plot(self):
        raise NotImplementedError
    
    def save(self):
        raise NotImplementedError


class HeatmapPlotter(Plotter):
    def __init__(self, matrix_path: str, target_path: str, plot_attr = {}) -> None:
        self.matrix_path = matrix_path
        self.target_path = target_path
        self.dict = {
            'xlabel': '', 
            'ylabel': '', 
            'fontsize': 22, 
            'figsize': (8, 5), 
            'cmap': "viridis"
            }
        # update selected attributes in dictionary
        for key in plot_attr:
            self.dict[key] = plot_attr[key]    


    def plot(self, show = True) -> None:
        # parse data from file
        with open(self.matrix_path, "r") as f:
            matrix = [float(line.split()[-1]) for line in f]
            matrix = np.array(matrix).reshape([int(np.sqrt(len(matrix))), -1])

        # create custom plot
        plt.ioff()
        fig = plt.figure(figsize = self.dict['figsize'])
        plt.imshow(matrix[:,3:-3], cmap = self.dict["cmap"])
        plt.colorbar()
        plt.xlabel(self.dict["xlabel"], fontsize = self.dict["fontsize"])
        plt.ylabel(self.dict["ylabel"], fontsize = self.dict["fontsize"])
        plt.tight_layout()
        if show:
            plt.show()


    def save(self) -> None:
        self.plot(show=False)

        plt.savefig(self.target_path)


    def update(self, new_figure_path='', new_target_path='', new_plot_attr={}) -> None:
        if new_figure_path: self.figure_path = new_figure_path
        if new_target_path: self.target_path = new_target_path
        if new_plot_attr:
            for key in new_plot_attr:
                self.dict[key] = new_plot_attr[key]


class ErrorPlotter(Plotter):
    def __init__(self, x_data: list, y_data: list, y_error: list, x_error:list = None, target_path: str = '', plot_attr = {}) -> None:
        self.x = x_data
        self.y = y_data
        self.x_err = x_error
        self.y_err = y_error
        self.target_path = target_path
        self.dict = {
            'title': '',
            'xlabel': '', 
            'ylabel': '', 
            'fontsize': 22, 
            'figsize': (8, 5), 
            'c': "black",
            'ecolor': "black",
            'capsize': 4,
            'fmt':'o',
            'grid_style': '',
            'dpi': 300,
            'ylim': (0, 1)
            }
        # update selected attributes in dictionary
        for key in plot_attr:
            self.dict[key] = plot_attr[key]    


    def plot(self, show=True) -> None:
        plt.figure(figsize=self.dict['figsize'])
        plt.errorbar(self.x, self.y, xerr=self.x_err, yerr=self.y_err, 
                     fmt=self.dict['fmt'], capsize=self.dict['capsize'], 
                     c=self.dict['c'], ecolor = self.dict['ecolor'])
        plt.title(self.dict['title'])
        plt.xlabel(self.dict['xlabel'], fontsize=self.dict['fontsize'])
        plt.ylabel(self.dict['ylabel'], fontsize=self.dict['fontsize'])
        plt.ylim(self.dict['ylim'][0], self.dict['ylim'][1])
        plt.tight_layout()

        if self.dict['grid_style']: 
            plt.grid(linestyle=self.dict['grid_style'])
        
        if show:
            plt.show()


    def save(self) -> None:
        self.plot(show=False)   

        plt.savefig(self.target_path, dpi = self.dict['dpi'])


    def update(self) -> None:
        # currently not needed
        pass