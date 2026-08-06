/**
 * spo2_pipeline.h
 *
 * SpO2 computation pipeline from PPG signals (IR + Red).
 *
 * Method: AC extraction via narrow-band DFT around the dominant heart-rate
 * frequency. Robust at low perfusion index (<1%).
 *
 * Requirements: C99, libm (-lm)
 * No dynamic allocation. All buffers are statically sized.
 */

#ifndef SPO2_PIPELINE_H
#define SPO2_PIPELINE_H

#include <stdint.h>

/* ─── Configuration ──────────────────────────────────────────────────── */

/* Sampling frequency of the IR/Red capture, in Hz. */
#define SPO2_FS 50.0f
/* Upper bound on the number of samples per capture (per channel) that the
 * statically-sized working buffers can hold. The actual sample count is
 * variable and always taken from the caller's raw capture length (the
 * n_samples argument to spo2_compute_R / spo2_compute); this macro only
 * caps how large that buffer can be, it is not the capture length itself. */
#define SPO2_MAX_N_SAMPLES 2000
/* Samples discarded at each edge of the capture, to crop out the bandpass
 * filter's startup transient (see bandpass_filter in spo2_pipeline.c). */
#define SPO2_MARGIN_SAMPLES 25

/* Lower edge of the cardiac frequency band, in Hz. */
#define SPO2_FC_MIN 0.5f
/* Upper edge of the cardiac frequency band, in Hz (~240 bpm). */
#define SPO2_FC_MAX 4.0f
/* Half-width of the narrow-band DFT window placed around the detected
 * heart-rate frequency, in Hz (i.e. the window spans f_c ± SPO2_DFT_BW). */
#define SPO2_DFT_BW 0.1f
/* Minimum acceptable SNR at the heart-rate peak, in dB; captures scoring
 * below this are rejected as too noisy to trust. */
#define SPO2_SNR_THRESHOLD 9.0f
/* Lower bound of the physiologically plausible range for the ratio R. */
#define SPO2_R_MIN 0.1f
/* Upper bound of the physiologically plausible range for the ratio R.
 * Capped at 1.0 because only ~0.3% of the calibration dataset (2 of 591
 * valid captures) falls above it; regressing past that point is
 * extrapolation from almost no data, not a supported measurement. */
#define SPO2_R_MAX 2.0f

/* SQI (Signal Quality Index), range 0-255. Combines three components,
 * each normalized to [0,1] and weighted by SQI_W_*, then scaled to 8 bits:
 *   - SNR: linear ramp between SQI_SNR_MIN and SQI_SNR_MAX (dB)
 *   - PI:  trapezoid, penalizes both weak perfusion and abnormally high PI
 *          (possible motion artifact / sensor decoupling)
 *   - Spectral purity: fraction of the cardiac-band energy concentrated at
 *     the peak (vs. spread out into broadband noise) */
/* Below this SNR (dB), the SNR component of the SQI is 0. */
#define SQI_SNR_MIN 0.0f
/* At or above this SNR (dB), the SNR component of the SQI is 1. */
#define SQI_SNR_MAX 20.0f

/* PI thresholds derived from the P5/P25/P90/P99 percentiles of the
 * calibration dataset (673 captures, inputs/rd.csv). The PI scale in this
 * pipeline (AC/DC via Goertzel) is ~100x smaller than typical clinical
 * literature values (0.5-5%), so do NOT reuse textbook/manual thresholds
 * here. Update together with SPO2_CAL_* after every calibration session. */
/* Below this PI (%), the PI component of the SQI is 0. */
#define SQI_PI_LOW_MIN 0.003f
/* At or above this PI (%), the PI component of the SQI reaches 1. */
#define SQI_PI_LOW_OPT 0.005f
/* Up to this PI (%), the PI component of the SQI stays at 1. */
#define SQI_PI_HIGH_OPT 0.020f
/* At or above this PI (%), the PI component of the SQI is 0. */
#define SQI_PI_HIGH_MAX 0.040f

/* Weight of the SNR component in the SQI blend. */
#define SQI_W_SNR 0.4f
/* Weight of the PI component in the SQI blend. */
#define SQI_W_PI 0.3f
/* Weight of the spectral-purity component in the SQI blend. */
#define SQI_W_PURITY 0.3f

/* Quadratic calibration coefficients (SpO2 = CAL_A + CAL_B * R + CAL_C * R^2).
 * Degree-2 polynomial regression over the calibration captures with
 * R <= SPO2_R_MAX (589 of 591 valid captures). Unconstrained regression over
 * the full dataset turns upward past R ~= 1.06 (not physiologically valid),
 * which is why SPO2_R_MAX excludes that region instead of forcing the curve
 * flat to cover it: within [SPO2_R_MIN, SPO2_R_MAX] this fit is already
 * monotonically decreasing on its own (vertex at R ~= 1.35, outside range).
 * RMSE = 4.31%
 * Update after every calibration session. */
#define SPO2_CAL_A 107.49f
#define SPO2_CAL_B -52.88f
#define SPO2_CAL_C 19.55f

/* Physiologically plausible bounds for the final SpO2 output (%); the
 * quadratic fit overshoots 100 near SPO2_R_MIN, so the output is clamped
 * to this range regardless of what the raw polynomial evaluates to. */
#define SPO2_OUT_MIN 70.0f
#define SPO2_OUT_MAX 100.0f

/* ─── Return codes ───────────────────────────────────────────────────── */

typedef enum
{
    SPO2_OK = 0,
    /* AC amplitude too low to extract a meaningful ratio (no detectable
     * pulsatile signal on the IR channel). */
    SPO2_ERR_LOW_AC = -1,
    /* SNR at the heart-rate peak below SPO2_SNR_THRESHOLD. */
    SPO2_ERR_LOW_SNR = -2,
    /* Computed R outside [SPO2_R_MIN, SPO2_R_MAX]. */
    SPO2_ERR_R_RANGE = -3,
    /* n_samples is too small to leave a usable window after cropping
     * SPO2_MARGIN_SAMPLES off each edge, or exceeds SPO2_MAX_N_SAMPLES. */
    SPO2_ERR_N_SAMPLES = -4,
} Spo2Status;

/* ─── Public API ─────────────────────────────────────────────────────── */

/**
 * Computes the ratio R and signal-quality metrics.
 *
 * @param ir        Array of n_samples IR samples (infrared channel)
 * @param red       Array of n_samples Red samples (red channel)
 * @param n_samples Number of samples per channel in this capture, as taken
 *                  from the raw signal (must be <= SPO2_MAX_N_SAMPLES and
 *                  leave a non-empty window after cropping
 *                  SPO2_MARGIN_SAMPLES off each edge, otherwise
 *                  SPO2_ERR_N_SAMPLES is returned)
 * @param out_R       [out] Ratio R = (AC_red/DC_red) / (AC_ir/DC_ir)
 * @param out_snr     [out] SNR at the dominant heart-rate frequency (dB)
 * @param out_pi      [out] IR perfusion index (%)
 * @param out_hr_bpm  [out] Dominant heart rate (beats per minute, integer)
 * @param out_sqi     [out] Signal Quality Index, 0-255 (255 = best quality)
 * @return            SPO2_OK, or an error code
 */
Spo2Status spo2_compute_R(const float *ir, const float *red, int n_samples,
                          float *out_R, float *out_snr, float *out_pi,
                          float *out_hr_bpm, uint8_t *out_sqi);

/**
 * Applies the quadratic calibration curve to obtain SpO2 (%).
 * SpO2 = CAL_A + CAL_B * R + CAL_C * R^2
 *
 * @param R     Valid ratio R (from spo2_compute_R)
 * @return      Estimated SpO2 (%)
 */
float spo2_predict(float R);

/**
 * Convenience function: compute_R + predict in a single call.
 * Returns SPO2_OK and fills *out_spo2, *out_hr_bpm and *out_sqi, or an
 * error code (in which case *out_spo2 is left untouched).
 */
Spo2Status spo2_compute(const float *ir, const float *red, int n_samples,
                        float *out_spo2, float *out_hr_bpm, uint8_t *out_sqi);

#endif /* SPO2_PIPELINE_H */
