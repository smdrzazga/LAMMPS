import numpy as np
import banana_lib as sz
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



plt.rcParams.update({'font.size': 18})
plt.rcParams.update({"text.usetex": False})


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
            'fontsize': 18, 
            'figsize': (6.4, 4.8), 
            'cmap': "viridis",
            'extent': (),
            'cbarlabel': ''
            }
        # update selected attributes in dictionary
        for key in plot_attr:
            self.dict[key] = plot_attr[key]    


    def plot(self, fig: plt.figure, ax: plt.axes, show = True) -> plt.figure:
        # parse data from file
        with open(self.matrix_path, "r") as f:
            matrix = [float(line.split()[-1]) for line in f]
            matrix = np.array(matrix).reshape([int(np.sqrt(len(matrix))), -1])

        # draw on axes and add colorbar on canvas eg. figure
        pos = ax.imshow(matrix[:,3:-3], cmap = self.dict["cmap"], extent=self.dict['extent'], aspect=70/45)
        fig.colorbar(pos, label=self.dict['cbarlabel'])
        ax.set_xlabel(self.dict["xlabel"], fontsize = self.dict["fontsize"])
        ax.set_ylabel(self.dict["ylabel"], fontsize = self.dict["fontsize"])
        plt.tight_layout()
        
        if show:
            ax.show()

        return ax


    def save(self) -> None:
        fig, ax = plt.subplots(figsize = self.dict["figsize"])
        self.plot(fig, ax, show=False)

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
            'label': '',
            'xlabel': '', 
            'ylabel': '', 
            'text': '',
            'textpos': (1, 1),
            'fontsize': 22, 
            'figsize': (8, 5), 
            'c': "black",
            'ecolor': "black",
            'ctick': 'black',
            'cylabel': 'black',
            'capsize': 4,
            'fmt':'o',
            'grid_style': '',
            'dpi': 300,
            'ylim': (0, 1)
            }
        # update selected attributes in dictionary
        for key in plot_attr:
            self.dict[key] = plot_attr[key]    


    def plot(self, ax: plt.axes, show=True) -> plt.figure:
        # fig = plt.figure(figsize=self.dict['figsize'])
        ax.errorbar(self.x, self.y, xerr=self.x_err, yerr=self.y_err, 
                     fmt=self.dict['fmt'], capsize=self.dict['capsize'], 
                     c=self.dict['c'], ecolor = self.dict['ecolor'], label=self.dict["label"])
        ax.set_title(self.dict['title'])
        ax.set_xlabel(self.dict['xlabel'], fontsize=self.dict['fontsize'])
        ax.set_ylabel(self.dict['ylabel'], fontsize=self.dict['fontsize'], c=self.dict['cylabel'])
        ax.tick_params('y', colors=self.dict['ctick'])
        ax.set_ylim(self.dict['ylim'][0], self.dict['ylim'][1])
        # ax.text(self.dict['textpos'][0], self.dict['textpos'][1], self.dict['text'])
        ax.text(self.dict['textpos'][0], self.dict['textpos'][1], self.dict['text'], bbox={'alpha': 0.2, 'color': '#D63230'})
        plt.tight_layout()

        if self.dict['grid_style']: 
            ax.grid(linestyle=self.dict['grid_style'])
        
        if show:
            ax.show()

        return ax


    def save(self) -> None:
        fig = plt.figure(figsize = self.dict["figsize"])
        self.plot(fig, show=False)   

        plt.savefig(self.target_path, dpi=self.dict['dpi'])


    def update(self) -> None:
        # currently not needed
        pass


class ScatterPlotter(Plotter):
    def __init__(self, matrix_path: str, target_path: str, plot_attr = {}) -> None:
        self.matrix_path = matrix_path
        self.target_path = target_path
        self.dict = {
            'xlabel': '', 
            'ylabel': '', 
            'fontsize': 18, 
            'c': 'black',
            's': 1/72,
            'figsize': (6, 6), 
            'xlim': (-0.7, 0.7),
            'ylim': (-0.7, 0.7),
            'ecolor': 'black',
            'label': ''
            }
        # update selected attributes in dictionary
        for key in plot_attr:
            self.dict[key] = plot_attr[key]    


    def plot(self, fig: plt.figure, ax: plt.axes, show = True, i: int = 1) -> plt.figure:
        # parse data from file
        with open(self.matrix_path, "r") as f:
            matrix = [np.array(line.split()[-3:], dtype=np.float64) for line in f]
            matrix = np.array(matrix).reshape([-1, 3])

        # draw on axes and add colorbar on canvas eg. figure
        x = matrix[:, 0] + 0.3*(i-1)
        y = matrix[:, 1] 
        pos = ax.scatter(x, y, c = self.dict["c"], s = self.dict['s'], edgecolors = self.dict['ecolor'])
        pos.set_label(self.dict['label'])
        
        # specify additional params
        ax.set_xlim(*self.dict["xlim"])
        ax.set_ylim(*self.dict["ylim"])
        ax.set_xlabel(self.dict["xlabel"], fontsize = self.dict["fontsize"])
        ax.set_ylabel(self.dict["ylabel"], fontsize = self.dict["fontsize"])
        plt.tight_layout()
        
        if show:
            plt.show()

        return ax


    def save(self) -> None:
        fig, ax = plt.subplots(figsize = self.dict["figsize"])
        self.plot(fig, ax, show=False)

        plt.savefig(self.target_path)


    def update(self, new_figure_path='', new_target_path='', new_plot_attr={}) -> None:
        if new_figure_path: self.figure_path = new_figure_path
        if new_target_path: self.target_path = new_target_path
        if new_plot_attr:
            for key in new_plot_attr:
                self.dict[key] = new_plot_attr[key]