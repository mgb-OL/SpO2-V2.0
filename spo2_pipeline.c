/**
 * spo2_pipeline.c
 *
 * Implementación del pipeline SpO2.
 * Compilar con: gcc -O2 -std=c99 -lm spo2_pipeline.c -o spo2
 */

#include "spo2_pipeline.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

/* ═══════════════════════════════════════════════════════════════════════
 * 1. FILTRO PASO BANDA  (Butterworth orden 4, Direct Form II transpuesta)
 *
 * Coeficientes SOS generados con:
 *   scipy.signal.butter(4, [0.5/25, 4.0/25], btype='band', output='sos')
 * fs=50 Hz, fl=0.5 Hz, fh=4.0 Hz  →  4 secciones biquad
 * ═══════════════════════════════════════════════════════════════════════ */

#define N_BIQUAD 4

/* Cada fila: { b0, b1, b2, a0(=1), a1, a2 } */
static const float SOS[N_BIQUAD][6] = {
    {1.3974753594e-03f, 2.7949507187e-03f, 1.3974753594e-03f,
     1.0000000000e+00f, -1.3776426634e+00f, 4.9806154978e-01f},
    {1.0000000000e+00f, 2.0000000000e+00f, 1.0000000000e+00f,
     1.0000000000e+00f, -1.5416687975e+00f, 7.4855141536e-01f},
    {1.0000000000e+00f, -2.0000000000e+00f, 1.0000000000e+00f,
     1.0000000000e+00f, -1.8683105055e+00f, 8.7412929192e-01f},
    {1.0000000000e+00f, -2.0000000000e+00f, 1.0000000000e+00f,
     1.0000000000e+00f, -1.9575073847e+00f, 9.6156538689e-01f},
};

/**
 * Aplica el filtro paso banda (causal, un solo paso).
 * Para uso en dispositivo: filtfilt no es necesario si se descarta
 * el transitorio inicial (ya cubierto por SPO2_MARGIN_SAMPLES).
 *
 * @param in    señal de entrada, longitud n
 * @param out   señal filtrada, longitud n (puede ser el mismo buffer)
 * @param n     número de muestras
 */
static void bandpass_filter(const float *in, float *out, int n)
{
    /* Estado de los delays para cada sección */
    float w[N_BIQUAD][2];
    memset(w, 0, sizeof(w));

    for (int i = 0; i < n; i++)
    {
        float x = in[i];
        for (int s = 0; s < N_BIQUAD; s++)
        {
            float b0 = SOS[s][0], b1 = SOS[s][1], b2 = SOS[s][2];
            float a1 = SOS[s][4], a2 = SOS[s][5];
            /* Direct Form II transpuesta */
            float y = b0 * x + w[s][0];
            w[s][0] = b1 * x - a1 * y + w[s][1];
            w[s][1] = b2 * x - a2 * y;
            x = y;
        }
        out[i] = x;
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * 2. DFT PUNTUAL  (Goertzel en banda cardíaca)
 *
 * En lugar de una FFT completa, calculamos la DFT solo en las frecuencias
 * de interés (0.5–4.0 Hz en pasos de 0.1 Hz). Esto es O(N·K) con K≈35,
 * eficiente en microcontrolador sin librería FFT externa.
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * Algoritmo de Goertzel: potencia espectral en frecuencia f.
 * Incluye ventana de Hanning para reducir spectral leakage.
 *
 * @param x     señal (con DC ya restado)
 * @param n     longitud
 * @param f     frecuencia objetivo (Hz)
 * @param fs    frecuencia de muestreo (Hz)
 * @return      potencia (unidades al cuadrado)
 */
static float goertzel_power(const float *x, int n, float f, float fs)
{
    float omega = 2.0f * (float)M_PI * f / fs;
    float coeff = 2.0f * cosf(omega);
    float s_prev = 0.0f, s_prev2 = 0.0f;

    for (int i = 0; i < n; i++)
    {
        /* Ventana de Hanning */
        float w = 0.5f * (1.0f - cosf(2.0f * (float)M_PI * i / (n - 1)));
        float s = w * x[i] + coeff * s_prev - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    }
    /* Potencia normalizada */
    float power = (s_prev2 * s_prev2 + s_prev * s_prev - coeff * s_prev * s_prev2) / (float)(n * n);
    return power;
}

/**
 * Encuentra la frecuencia cardíaca dominante (Hz) en [FC_MIN, FC_MAX].
 * Resolución: 0.1 Hz.
 */
static float find_dominant_freq(const float *ac_signal, int n, float fs)
{
    float best_power = -1.0f;
    float best_freq = SPO2_FC_MIN;

    for (float f = SPO2_FC_MIN; f <= SPO2_FC_MAX + 0.01f; f += 0.1f)
    {
        float p = goertzel_power(ac_signal, n, f, fs);
        if (p > best_power)
        {
            best_power = p;
            best_freq = f;
        }
    }
    return best_freq;
}

/**
 * Amplitud RMS en la banda [f_c - BW, f_c + BW].
 * Suma la potencia en los bins de la ventana y devuelve sqrt.
 */
static float dft_amplitude(const float *signal_no_dc, int n, float fs,
                           float f_c, float bw)
{
    float power_sum = 0.0f;
    int count = 0;

    for (float f = f_c - bw; f <= f_c + bw + 0.01f; f += 0.1f)
    {
        if (f < 0.01f)
            continue;
        power_sum += goertzel_power(signal_no_dc, n, f, fs);
        count++;
    }
    return (count > 0) ? sqrtf(power_sum) : 0.0f;
}

/**
 * SNR (dB) del pico cardíaco respecto al ruido de fondo en la banda.
 */
static float compute_snr(const float *ac_signal, int n, float fs,
                         float f_c, float bw)
{
    float peak_power = 0.0f;
    float total_power = 0.0f;
    int total_bins = 0;

    for (float f = SPO2_FC_MIN; f <= SPO2_FC_MAX + 0.01f; f += 0.1f)
    {
        float p = goertzel_power(ac_signal, n, f, fs);
        total_power += p;
        total_bins++;
        if (f >= f_c - bw && f <= f_c + bw)
        {
            if (p > peak_power)
                peak_power = p;
        }
    }

    float noise_mean = (total_bins > 0)
                           ? total_power / (float)total_bins
                           : 1e-12f;

    if (noise_mean < 1e-30f)
        return 0.0f;
    return 10.0f * log10f(peak_power / noise_mean);
}

/* ═══════════════════════════════════════════════════════════════════════
 * 3. CÁLCULO DEL RATIO R
 * ═══════════════════════════════════════════════════════════════════════ */

Spo2Status spo2_compute_R(const float *ir, const float *red,
                          float *out_R, float *out_snr, float *out_pi)
{
    /* Buffers estáticos: sin malloc */
    static float ir_work[SPO2_N_WORK];
    static float red_work[SPO2_N_WORK];
    static float ir_ac[SPO2_N_WORK];
    static float ir_nodC[SPO2_N_WORK];
    static float red_nodC[SPO2_N_WORK];

    const int m = SPO2_MARGIN_SAMPLES;
    const int n = SPO2_N_WORK;

    /* Copiar región de trabajo (descartamos bordes) */
    memcpy(ir_work, ir + m, n * sizeof(float));
    memcpy(red_work, red + m, n * sizeof(float));

    /* DC = media aritmética */
    float dc_ir = 0.0f, dc_red = 0.0f;
    for (int i = 0; i < n; i++)
    {
        dc_ir += ir_work[i];
        dc_red += red_work[i];
    }
    dc_ir /= n;
    dc_red /= n;

    /* Señal sin DC para DFT */
    for (int i = 0; i < n; i++)
    {
        ir_nodC[i] = ir_work[i] - dc_ir;
        red_nodC[i] = red_work[i] - dc_red;
    }

    /* Filtrado paso banda para detectar FC dominante */
    bandpass_filter(ir_nodC, ir_ac, n);

    /* Frecuencia cardíaca dominante */
    float f_c = find_dominant_freq(ir_ac, n, SPO2_FS);

    /* Amplitud AC via DFT en banda estrecha */
    float amp_ir = dft_amplitude(ir_nodC, n, SPO2_FS, f_c, SPO2_DFT_BW);
    float amp_red = dft_amplitude(red_nodC, n, SPO2_FS, f_c, SPO2_DFT_BW);

    if (amp_ir < 1e-6f)
    {
        *out_R = 0.0f;
        *out_snr = 0.0f;
        *out_pi = 0.0f;
        return SPO2_ERR_LOW_AC;
    }

    float R = (amp_red / dc_red) / (amp_ir / dc_ir);
    float snr = compute_snr(ir_ac, n, SPO2_FS, f_c, SPO2_DFT_BW);
    float pi = (amp_ir / dc_ir) * 100.0f;

    *out_R = R;
    *out_snr = snr;
    *out_pi = pi;

    /* Validación de calidad */
    if (snr < SPO2_SNR_THRESHOLD)
        return SPO2_ERR_LOW_SNR;
    if (R < SPO2_R_MIN || R > SPO2_R_MAX)
        return SPO2_ERR_R_RANGE;

    return SPO2_OK;
}

/* ═══════════════════════════════════════════════════════════════════════
 * 4. PREDICCIÓN DE SpO2
 * ═══════════════════════════════════════════════════════════════════════ */

float spo2_predict(float R)
{
    return SPO2_CAL_A + SPO2_CAL_B * R + SPO2_CAL_C * R * R;
}

Spo2Status spo2_compute(const float *ir, const float *red, float *out_spo2)
{
    float R, snr, pi;
    Spo2Status status = spo2_compute_R(ir, red, &R, &snr, &pi);
    if (status == SPO2_OK)
    {
        *out_spo2 = spo2_predict(R);
    }
    return status;
}
