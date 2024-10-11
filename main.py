from common.PROCESSING import *
from common.TIMERS import *
from common.IO import *
from config import GlobalParameters



#TODO: move timer to analyzer

if __name__ == '__main__':
    params = GlobalParameters()
    analyzer = DirectorFullAnalyzer(params)

    timer_analysis = Timer()
    timer_printer = Timer(message_start='Here comes the heatmap!')

    timer_analysis.print_duration(
        analyzer.get_results_in_parallel()
    )

    screen = analyzer.screen
    timer_printer.print_duration(
        ScreenPrinter(screen).print_screen(params['TARGET_FILE'])
    )
