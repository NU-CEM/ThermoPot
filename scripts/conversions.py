"""Script converts from older version of phonopy output to newer version.
Both the units and the column order have changed.

Older version has data organised as: T(K) F(kJ/mol)  S(J/K/mol) Cv(J/Kmol)  U(kJ/mol).
Newer version has data organised as: T(K) F(eV/cell) U(eV/cell) Cv(kB/cell) -TS(eV/cell)."""


def convert_phonopy_filetypes(file_in, file_out):

    df_in = pd.read_csv(
        file_in, comment="#", delim_whitespace=True, names=["T", "F", "U", "Cv", "-TS"]
    )
    df_out = df_in
    df_out["F"] = df_in["F"] * thermopot.eV2kJmol
    df_out["U"] = df_in["U"] * thermopot.eV2kJmol
    df_out["Cv"] = df_in["Cv"] * thermopot.kB2JKmol
    df_out["-TS"] = df_in["-TS"] * thermopot.kB2JKmol / (-1 * df_in["T"])
    df_out = df_out.rename({"-TS": "S"}, axis=1)
    df_out = df_out[["T", "F", "S", "Cv", "U"]]
    df_out = df_out.fillna(0)
    df_out.to_csv(file_out, sep=" ", index=False, header=False)

    return df_out
