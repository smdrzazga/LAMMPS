from common.SCREEN import *
from common.ANALYZERS import *


N = 150
N_SLICES = 75
SMECTIC_PERIODS = 4

if __name__ == "__main__":
    analyzer = SmecticAnalyzer(N_SLICES, SMECTIC_PERIODS=SMECTIC_PERIODS, N_PIX=N)
    analyzer.print_smectic_params_for_all_files("centers_screen_bulk_d_0.31_*")

    analyzer = CorrelationLengthAnalyzer()
    analyzer.print_correlation_params_for_all_files("smectic_params_0.31*")

    analyzer = EllipsisAnalyzer(N_SLICES, N_PIX=N)
    analyzer.print_semiaxes_for_all_files("polarization_screen_full_d_0.31*")
    