from common.PROCESSING import *
from common.TIMERS import *
from common.IO import *
from config import GlobalParameters



#TODO: move timer to analyzer

if __name__ == '__main__':
    input_file = 'C:/Users/Szymek/Desktop/middle_snapshot_4000000.lammpstrj'
    params = GlobalParameters(INPUT_FILE = input_file)
    analyzer = DirectorFullAnalyzer(params)

    # timer_analysis = Timer()
    # timer_printer = Timer(message_start='Here comes the heatmap!')

    # analyzer.get_results_in_parallel()
    analyzer.get_results_in_parallel()
    
    # screen = analyzer.screen
    # ScreenPrinter(screen).print_screen(params.file_params['OUTPUT_FILE'])

