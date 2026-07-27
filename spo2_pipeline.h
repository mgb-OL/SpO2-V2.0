/**
 * spo2_pipeline.h
 *
 * Pipeline de cálculo de SpO2 a partir de señales PPG (IR + Red).
 *
 * Método: extracción de AC vía DFT en banda estrecha alrededor de la
 * frecuencia cardíaca dominante. Robusto con perfusion index bajo (<1%).
 *
 * Requisitos: C99, libm (-lm)
 * Sin allocación dinámica. Todos los buffers son de tamaño estático.
 */

#ifndef SPO2_PIPELINE_H
#define SPO2_PIPELINE_H

#include <stdint.h>

/* ─── Configuración ─────────────────────────────────────────────────── */

#define SPO2_FS 50.0f          /* Hz — frecuencia de muestreo        */
#define SPO2_N_SAMPLES 600     /* muestras por captura               */
#define SPO2_MARGIN_SAMPLES 25 /* muestras a descartar en cada borde */
#define SPO2_N_WORK 550        /* N_SAMPLES - 2*MARGIN               */

#define SPO2_FC_MIN 0.5f        /* Hz — límite inferior banda cardíaca */
#define SPO2_FC_MAX 4.0f        /* Hz — límite superior (~240 bpm)     */
#define SPO2_DFT_BW 0.1f        /* Hz — semiancho ventana DFT          */
#define SPO2_SNR_THRESHOLD 9.0f /* dB — capturas por debajo se rechazan */
#define SPO2_R_MIN 0.1f         /* límite fisiológico inferior de R     */
#define SPO2_R_MAX 2.0f         /* límite fisiológico superior de R     */

/* Coeficientes de calibración cuadrática (SpO2 = CAL_A + CAL_B * R + CAL_C * R²)
 * Obtenidos por regresión polinomial de grado 2 sobre dataset de calibración.
 * RMSE = 4.27%  (vs 4.49% del modelo lineal)
 * Actualizar tras cada sesión de calibración. */
#define SPO2_CAL_A 108.84f
#define SPO2_CAL_B -60.25f
#define SPO2_CAL_C 28.34f

/* ─── Códigos de retorno ─────────────────────────────────────────────── */

typedef enum
{
    SPO2_OK = 0,
    SPO2_ERR_LOW_AC = -1,  /* amplitud AC demasiado baja             */
    SPO2_ERR_LOW_SNR = -2, /* SNR por debajo del umbral              */
    SPO2_ERR_R_RANGE = -3, /* R fuera del rango fisiológico          */
} Spo2Status;

/* ─── API pública ────────────────────────────────────────────────────── */

/**
 * Calcula el ratio R y métricas de calidad de señal.
 *
 * @param ir        Array de SPO2_N_SAMPLES muestras IR  (canal infrarrojo)
 * @param red       Array de SPO2_N_SAMPLES muestras Red (canal rojo)
 * @param out_R     [out] Ratio R = (AC_red/DC_red) / (AC_ir/DC_ir)
 * @param out_snr   [out] SNR en la frecuencia cardíaca dominante (dB)
 * @param out_pi    [out] Perfusion index IR (%)
 * @return          SPO2_OK o código de error
 */
Spo2Status spo2_compute_R(const float *ir, const float *red,
                          float *out_R, float *out_snr, float *out_pi);

/**
 * Aplica la curva de calibración cuadrática para obtener SpO2 (%).
 * SpO2 = CAL_A + CAL_B * R + CAL_C * R²
 *
 * @param R     Ratio R válido (de spo2_compute_R)
 * @return      SpO2 estimado (%)
 */
float spo2_predict(float R);

/**
 * Función de conveniencia: compute_R + predict en un solo paso.
 * Devuelve SPO2_OK y rellena *out_spo2, o un código de error.
 */
Spo2Status spo2_compute(const float *ir, const float *red, float *out_spo2);

#endif /* SPO2_PIPELINE_H */
