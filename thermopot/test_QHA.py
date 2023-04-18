from thermopot import materials, calculations, reactions
import numpy as np

BaZrS3_calc = calculations.QHACalculation(
    filepath="../BaZrS3/raw_aims_files/ternary/BaZrS3_Pnma/hse06/aims.out",
    ev_filepath="../BaZrS3/qha_output/BaZrS3_Pnma/e-v.dat",
)

Ba2ZrS4_calc = calculations.QHACalculation(
    filepath="../BaZrS3/raw_aims_files/ternary/Ba2ZrS4_I4_mmm/hse06/aims.out",
    ev_filepath="../BaZrS3/qha_output/Ba2ZrS4_I4_mmm/e-v.dat",
)

ZrS2_calc = calculations.QHACalculation(
    filepath="../BaZrS3/raw_aims_files/binary/ZrS2_P-3m1/hse06/aims.out",
    ev_filepath="../BaZrS3/qha_output/ZrS2_P-3m1/e-v.dat",
)

BaZrS3 = materials.Solid(
    "BaZrS3",
    {"Ba": 1, "Zr": 1, "S": 3},
    phonon_filepath="../BaZrS3/phonopy_output/BaZrS3_Pnma.dat",
    qha_calculation=BaZrS3_calc,
)

Ba2ZrS4 = materials.Solid(
    "Ba2ZrS4",
    {"Ba": 2, "Zr": 1, "S": 4},
    phonon_filepath="../BaZrS3/phonopy_output/Ba2ZrS4_I4_mmm.dat",
    qha_calculation=Ba2ZrS4_calc,
)

ZrS2 = materials.Solid(
    "ZrS2",
    {"Zr": 1, "S": 2},
    phonon_filepath="../BaZrS3/phonopy_output/ZrS2_P-3m1.dat",
    qha_calculation=ZrS2_calc,
)

# print(
#    BaZrS3.mu(
#        P=1,
#        T=300,
#        xc="hse06",
#        units="eV",
#    )
# )
T = np.linspace(100, 1300, 10)  # K
P = np.array(np.logspace(-3, 6, 100), ndmin=2).transpose()  # Pa
print(
    Ba2ZrS4.mu(
        P=P,
        T=T,
        xc="hse06",
        units="eV",
    )
)

# print(T.shape)
# reaction_214 = reactions.Reaction({BaZrS3:2},{Ba2ZrS4:1,ZrS2:1}, temperature=T, pressure=P, fu=2)
# print(vars(reaction_214))
# GFE_214 = reaction_214.Dmu(units="kJ",xc='hse06')
# GFE_214.plot_TvsP(scale_range=[17,15])
# GFE_214.potential[0][-1]
