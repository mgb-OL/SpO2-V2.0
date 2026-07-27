"""
Convierte inputs/rd.csv (celdas con arrays IR/Red embebidos como texto tipo
diccionario Python) en un CSV plano de columnas numericas:

    IR_0, ..., IR_599, Red_0, ..., Red_599, SpO2_ref

listo para que main.c (que enlaza con spo2_pipeline.c) lo lea directamente
con fscanf/strtok, sin necesidad de parsear texto en C.

Solo se incluyen capturas cuyos arrays IR y Red tengan exactamente
N_SAMPLES muestras (deben coincidir con SPO2_N_SAMPLES en spo2_pipeline.h).
"""
import sys

import pandas as pd

from spo2_pipeline import _extract_arrays

N_SAMPLES = 600  # debe coincidir con SPO2_N_SAMPLES en spo2_pipeline.h


def convert(csv_in: str, csv_out: str) -> int:
    df_raw = pd.read_csv(csv_in)
    rows = []

    for i in range(len(df_raw)):
        for j in range(len(df_raw.columns)):
            cell = df_raw.iloc[i, j]
            if not isinstance(cell, str):
                continue
            ir, red, ref = _extract_arrays(cell)
            if ir is None or red is None or ref is None:
                continue
            if len(ir) != N_SAMPLES or len(red) != N_SAMPLES:
                continue
            rows.append(list(ir) + list(red) + [ref])

    columns = ([f"IR_{k}" for k in range(N_SAMPLES)]
               + [f"Red_{k}" for k in range(N_SAMPLES)]
               + ["SpO2_ref"])
    pd.DataFrame(rows, columns=columns).to_csv(csv_out, index=False)
    return len(rows)


if __name__ == "__main__":
    csv_in = sys.argv[1] if len(sys.argv) > 1 else "inputs/rd.csv"
    csv_out = sys.argv[2] if len(sys.argv) > 2 else "inputs/spo2_captures.csv"

    n = convert(csv_in, csv_out)
    print(f"{n} capturas de {N_SAMPLES} muestras escritas en {csv_out}")
