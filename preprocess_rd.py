"""
Convierte inputs/rd.csv (celdas con arrays IR/Red embebidos como texto tipo
diccionario Python) en un CSV plano de columnas numericas, una fila por
captura:

    n_samples, IR_0, ..., IR_{n_samples-1}, Red_0, ..., Red_{n_samples-1},
    SpO2_ref, rd_row, rd_column

listo para que main.c (que enlaza con spo2_pipeline.c) lo lea directamente
con strtok, sin necesidad de parsear texto en C. n_samples es la longitud
real de cada captura (tomada del propio array IR/Red) y puede variar de una
fila a otra: no se asume ningun tamano fijo, cada fila lleva el suyo por
delante para que main.c sepa cuantos valores leer.

rd_row y rd_column identifican la celda de origen en rd.csv (fila y nombre
de columna, p.ej. "S3"), para poder cruzar los resultados de main.c con el
dataset original.

Solo se incluyen capturas cuyos arrays IR y Red tengan la misma longitud
entre si (y al menos una muestra); no se exige que coincidan con ningun
N_SAMPLES fijo. main.c descarta las que no le quepan en su buffer estatico
(ver SPO2_MAX_N_SAMPLES en spo2_pipeline.h).
"""
import csv
import sys

import pandas as pd

from spo2_pipeline import _extract_arrays


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
            if len(ir) != len(red) or len(ir) == 0:
                continue
            rows.append([len(ir)] + list(ir) + list(red) + [ref, i, df_raw.columns[j]])

    with open(csv_out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["n_samples", "IR...", "Red...", "SpO2_ref", "rd_row", "rd_column"])
        writer.writerows(rows)
    return len(rows)


if __name__ == "__main__":
    csv_in = sys.argv[1] if len(sys.argv) > 1 else "inputs/rd.csv"
    csv_out = sys.argv[2] if len(sys.argv) > 2 else "inputs/spo2_captures.csv"

    n = convert(csv_in, csv_out)
    print(f"{n} capturas escritas en {csv_out}")
