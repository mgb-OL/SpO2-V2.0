/**
 * main.c
 *
 * Programa de prueba para spo2_pipeline.c: lee capturas IR/Red desde un CSV
 * plano (generado por preprocess_rd.py a partir de inputs/rd.csv) y ejecuta
 * el pipeline de calculo de SpO2 sobre cada una, comparando contra la
 * referencia y reportando metricas de error.
 *
 * Ademas escribe un CSV de resultados (row,column,R,SNR,PI,valid,SpO2_ref,
 * SpO2_C) con la salida real del algoritmo en C para cada medicion, para
 * poder cruzarlo con el dataset original inputs/rd.csv.
 *
 * Archivo autocontenido (incluye spo2_pipeline.c directamente): basta con
 * compilar/ejecutar este archivo solo, sin enlazar nada mas.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Se incluye la implementacion (no solo el header) para que este archivo
 * sea autocontenido: los botones "Run"/"Debug" de los runners de VS Code
 * compilan unicamente el archivo activo, sin enlazar otros .c del proyecto. */
#include "spo2_pipeline.c"

/* Una fila tiene 2*SPO2_N_SAMPLES valores + SpO2_ref + rd_row + rd_column,
 * cada uno hasta ~12 caracteres (incluida la coma); se deja margen amplio. */
#define MAX_LINE (SPO2_N_SAMPLES * 2 * 16 + 64)
#define MAX_COLUMN_LEN 32

static int parse_row(char *line, float *ir, float *red, float *ref,
                      int *rd_row, char *rd_column)
{
    char *tok = strtok(line, ",\r\n");

    for (int i = 0; i < SPO2_N_SAMPLES; i++)
    {
        if (!tok)
            return 0;
        ir[i] = strtof(tok, NULL);
        tok = strtok(NULL, ",\r\n");
    }
    for (int i = 0; i < SPO2_N_SAMPLES; i++)
    {
        if (!tok)
            return 0;
        red[i] = strtof(tok, NULL);
        tok = strtok(NULL, ",\r\n");
    }
    if (!tok)
        return 0;
    *ref = strtof(tok, NULL);

    tok = strtok(NULL, ",\r\n");
    if (!tok)
        return 0;
    *rd_row = atoi(tok);

    tok = strtok(NULL, ",\r\n");
    if (!tok)
        return 0;
    strncpy(rd_column, tok, MAX_COLUMN_LEN - 1);
    rd_column[MAX_COLUMN_LEN - 1] = '\0';

    return 1;
}

int main(int argc, char **argv)
{
    const char *path = (argc > 1) ? argv[1] : "inputs/spo2_captures.csv";
    const char *out_path = (argc > 2) ? argv[2] : "outputs/spo2_results.csv";

    FILE *f = fopen(path, "r");
    if (!f)
    {
        fprintf(stderr, "No se pudo abrir '%s'\n", path);
        return 1;
    }

    FILE *fout = fopen(out_path, "w");
    if (!fout)
    {
        fprintf(stderr, "No se pudo crear '%s'\n", out_path);
        fclose(f);
        return 1;
    }
    fprintf(fout, "row,column,R,SNR,PI,valid,SpO2_ref,SpO2_C\n");

    static char line[MAX_LINE];
    static float ir[SPO2_N_SAMPLES];
    static float red[SPO2_N_SAMPLES];
    static char rd_column[MAX_COLUMN_LEN];

    if (!fgets(line, sizeof(line), f)) /* saltar cabecera */
    {
        fprintf(stderr, "Archivo vacio: '%s'\n", path);
        fclose(f);
        fclose(fout);
        return 1;
    }

    int n_total = 0, n_ok = 0;
    int n_low_ac = 0, n_low_snr = 0, n_r_range = 0;
    double sum_sq_err = 0.0, sum_abs_err = 0.0;

    /* Acumuladores de RMSE por nivel de SpO2_ref (niveles enteros 0-100). */
    #define N_LEVELS 101
    static int level_n[N_LEVELS];
    static double level_sum_sq[N_LEVELS];

    while (fgets(line, sizeof(line), f))
    {
        float ref;
        int rd_row;
        if (!parse_row(line, ir, red, &ref, &rd_row, rd_column))
            continue;

        n_total++;

        float R, snr, pi;
        Spo2Status status = spo2_compute_R(ir, red, &R, &snr, &pi);
        int valid = (status == SPO2_OK);
        float spo2 = roundf(spo2_predict(R));

        fprintf(fout, "%d,%s,%.9g,%.9g,%.9g,%s,%.9g,%.9g\n",
                rd_row, rd_column, (double)R, (double)snr, (double)pi,
                valid ? "True" : "False", (double)ref, (double)spo2);

        switch (status)
        {
        case SPO2_OK:
        {
            n_ok++;
            double err = (double)spo2 - (double)ref;
            sum_sq_err += err * err;
            sum_abs_err += fabs(err);

            int level = (int)lround((double)ref);
            if (level >= 0 && level < N_LEVELS)
            {
                level_n[level]++;
                level_sum_sq[level] += err * err;
            }
            break;
        }
        case SPO2_ERR_LOW_AC:
            n_low_ac++;
            break;
        case SPO2_ERR_LOW_SNR:
            n_low_snr++;
            break;
        case SPO2_ERR_R_RANGE:
            n_r_range++;
            break;
        }
    }
    fclose(f);
    fclose(fout);

    double pct_ok = (n_total > 0) ? (100.0 * n_ok / n_total) : 0.0;

    printf("\nResultados sobre %d capturas:\n", n_total);
    if (n_ok > 0)
    {
        printf("\nRMSE por nivel de SpO2:\n");
        for (int level = 0; level < N_LEVELS; level++)
        {
            if (level_n[level] == 0)
                continue;
            double rmse_level = sqrt(level_sum_sq[level] / level_n[level]);
            printf(" %3d%%  n=%3d  RMSE = %.2f %%\n", level, level_n[level], rmse_level);
        }

        double rmse = sqrt(sum_sq_err / n_ok);
        double mae = sum_abs_err / n_ok;
        printf("\nRMSE = %.2f %%\n", rmse);
        printf("MAE  = %.2f %%\n\n", mae);
    }

    printf("Capturas procesadas      = %d muestras\n", n_total);
    printf("Capturas validas         = %d (%.1f %%) \n", n_ok, pct_ok);
    printf(" Rechazadas por AC bajo          = %d\n", n_low_ac);
    printf(" Rechazadas por SNR bajo         = %d\n", n_low_snr);
    printf(" Rechazadas por R fuera de rango = %d\n\n", n_r_range);

    printf("Resultados escritos en '%s'\n\n", out_path);

    return 0;
}
