import os
from common.SCREEN import *
from common.ANALYZERS import *
import glob


N = 150
N_SLICES = 75
SMECTIC_PERIODS = 4

if __name__ == "__main__":
    # analyzer = SmecticAnalyzer(N_SLICES, SMECTIC_PERIODS=SMECTIC_PERIODS, N_PIX=N)
    # analyzer.print_smectic_params_for_all_files()

    # analyzer = CorrelationLengthAnalyzer()
    # analyzer.print__correlation_params_for_all_files()

    analyzer = EllipsisAnalyzer(N_SLICES, N_PIX=N)
    analyzer.print_semiaxes_for_all_files()