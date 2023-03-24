import calculations, materials, reactions
import numpy as np

data_path = "/Users/w21013885/QHA_code/ThermoPot/BaZrS3/"

BaZrS3_calc = calculations.AimsCalculation(
    data_path + "raw_aims_files/ternary/BaZrS3_Pnma/hse06/aims.out"
)
Ba2Zr1S4_calc = calculations.AimsCalculation(
    data_path + "raw_aims_files/ternary/Ba2ZrS4_I4_mmm/hse06/aims.out"
)
Ba3Zr2S7_calc = calculations.AimsCalculation(
    data_path + "raw_aims_files/ternary/Ba3Zr2S7_I4_mmm/hse06/aims.out"
)
Ba4Zr3S10_calc = calculations.AimsCalculation(
    data_path + "raw_aims_files/ternary/Ba4Zr3S10_I4_mmm/hse06/aims.out"
)
ZrS2_calc = calculations.AimsCalculation(
    data_path + "raw_aims_files/binary/ZrS2_P-3m1/hse06/aims.out"
)

BaZrS3 = materials.Solid(
    "BaZrS3",
    {"Ba": 1, "Zr": 1, "S": 3},
    data_path + "phonopy_output/BaZrS3_Pnma.dat",
    calculation=BaZrS3_calc,
)
Ba2Zr1S4 = materials.Solid(
    "Ba2Zr1S4",
    {"Ba": 2, "Zr": 1, "S": 4},
    data_path + "phonopy_output/Ba2ZrS4_I4_mmm.dat",
    calculation=Ba2Zr1S4_calc,
)
Ba3Zr2S7 = materials.Solid(
    "Ba3Zr2S7",
    {"Ba": 3, "Zr": 2, "S": 7},
    data_path + "phonopy_output/Ba3Zr2S7_I4_mmm.dat",
    calculation=Ba3Zr2S7_calc,
)
Ba4Zr3S10 = materials.Solid(
    "Ba4Zr3S10",
    {"Ba": 4, "Zr": 3, "S": 10},
    data_path + "phonopy_output/Ba4Zr3S10_I4_mmm.dat",
    calculation=Ba4Zr3S10_calc,
)
ZrS2 = materials.Solid(
    "ZrS2",
    {"Zr": 1, "S": 2},
    data_path + "/phonopy_output/ZrS2_P-3m1.dat",
    calculation=ZrS2_calc,
)

print(vars(BaZrS3))
T = np.linspace(100, 1300, 100)  # K
P = np.array(np.logspace(-3, 6, 100), ndmin=2).transpose()  # Pa

reaction_214 = reactions.Reaction(
    {BaZrS3: 2}, {Ba2Zr1S4: 1, ZrS2: 1}, temperature=T, pressure=P, fu=2
)
# print(vars(reaction_214))
