from common.SCREEN import *


class HorizontalSlice(Screen):
    def __init__(self, screen: Screen, borders: list[list]) -> None:
        size = np.array([borders[i][1] - borders[i][0] for i in range(2)])
        super().__init__(size, screen._get_pixel_type())
        self.screen = screen.get_slice(borders)
        
    def slice_director(self) -> list:
        if not isinstance(self.component, DirectorPixel):
            raise Exception("Pixel has to be Director Pixel")

        return self.avg_director()


class VerticalSlice(Screen):
    def __init__(self, screen: Screen, borders: list[list]) -> None:
        size = np.array([borders[i][1] - borders[i][0] for i in range(2)])
        super().__init__(size, screen._get_pixel_type())
        self.screen = screen.get_slice(borders)
        
    def slice_director(self) -> list:
        if not isinstance(self.component, DirectorPixel):
            raise Exception("Pixel has to be Director Pixel")

        return self.avg_director()


class SmecticParameter():
    def __init__(self, periods: int) -> None:
        self.parameter = 0 + 0j
        self.count = 0
        self.periods = periods

    def __repr__(self) -> str:
       return f"{np.abs(self.parameter):.3f} | {self.count}"

    # only for vertical slices in y-z plane!
    def add_atom(self, center_coords: list, box: SimulationBox) -> None:
        self.parameter += np.exp(self.periods * 2*np.pi*1j * center_coords[-1] / box.z)
        self.count += 1

    def normalize(self) -> None:
        if self.count != 0:
            self.parameter /= self.count
    
    def read_screen(self, screen: Screen, start: int, end: int) -> None:
        HARD_CODED_LIMIT = 0

        for x in range(start, end):
            for y in range(HARD_CODED_LIMIT, screen.y):
                self.parameter += np.exp(self.periods * 2*np.pi*1j * y / screen.y) * screen.screen[y][x].colour()
                self.count += screen.screen[y][x].colour()
