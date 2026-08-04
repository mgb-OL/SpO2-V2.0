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

/* SQI (Signal Quality Index), rango 0-255. Combina tres componentes
 * normalizados a [0,1] con pesos SQI_W_*, escalado luego a 8 bits:
 *   - SNR: rampa lineal entre SQI_SNR_MIN y SQI_SNR_MAX (dB)
 *   - PI:  trapezoide, penaliza perfusión débil y también PI anómalamente
 *          alto (posible artefacto de movimiento / desacople del sensor)
 *   - Pureza espectral: fracción de la energía en banda cardíaca
 *     concentrada en el pico (vs. dispersa en ruido de banda ancha) */
#define SQI_SNR_MIN 0.0f  /* dB — por debajo, componente SNR = 0     */
#define SQI_SNR_MAX 20.0f /* dB — por encima, componente SNR = 1     */

/* Umbrales de PI derivados de los percentiles P5/P25/P90/P99 del dataset
 * de calibración (673 capturas, inputs/rd.csv). La escala de PI en este
 * pipeline (AC/DC vía Goertzel) es ~100x menor que los valores típicos de
 * la literatura clínica (0.5-5%), así que NO usar umbrales de manual/paper.
 * Actualizar junto con SPO2_CAL_* tras cada sesión de calibración. */
#define SQI_PI_LOW_MIN 0.003f  /* %  — por debajo, componente PI = 0      */
#define SQI_PI_LOW_OPT 0.005f  /* %  — a partir de aquí, componente PI = 1 */
#define SQI_PI_HIGH_OPT 0.020f /* %  — hasta aquí, componente PI = 1      */
#define SQI_PI_HIGH_MAX 0.040f /* %  — por encima, componente PI = 0      */

#define SQI_W_SNR 0.4f    /* peso del componente SNR                 */
#define SQI_W_PI 0.3f     /* peso del componente PI                  */
#define SQI_W_PURITY 0.3f /* peso del componente de pureza espectral */

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
 * @param out_R       [out] Ratio R = (AC_red/DC_red) / (AC_ir/DC_ir)
 * @param out_snr     [out] SNR en la frecuencia cardíaca dominante (dB)
 * @param out_pi      [out] Perfusion index IR (%)
 * @param out_hr_bpm  [out] Frecuencia cardíaca dominante (pulsaciones por minuto, entero)
 * @param out_sqi     [out] Signal Quality Index, 0-255 (255 = mejor calidad)
 * @return            SPO2_OK o código de error
 */
Spo2Status spo2_compute_R(const float *ir, const float *red,
                          float *out_R, float *out_snr, float *out_pi,
                          float *out_hr_bpm, uint8_t *out_sqi);

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
 * Devuelve SPO2_OK y rellena *out_spo2, *out_hr_bpm y *out_sqi, o un
 * código de error.
 */
Spo2Status spo2_compute(const float *ir, const float *red, float *out_spo2,
                        float *out_hr_bpm, uint8_t *out_sqi);

#endif /* SPO2_PIPELINE_H */
