/**
 * spo2_pipeline.c
 *
 * Implementation of the SpO2 pipeline.
 * Compile with: gcc -O2 -std=c99 -lm spo2_pipeline.c -o spo2
 */

#include "spo2_pipeline.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

/* ═══════════════════════════════════════════════════════════════════════
 * 1. BANDPASS FILTER  (4th-order Butterworth, transposed Direct Form II)
 *
 * SOS coefficients generated with:
 *   scipy.signal.butter(4, [0.5/25, 4.0/25], btype='band', output='sos')
 * fs=50 Hz, fl=0.5 Hz, fh=4.0 Hz  ->  4 biquad sections
 * ═══════════════════════════════════════════════════════════════════════ */

#define N_BIQUAD 4

/* Each row: { b0, b1, b2, a0(=1), a1, a2 } */
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
 * Applies the bandpass filter (causal, single pass).
 * For on-device use: filtfilt is not needed as long as the initial
 * transient is discarded (already covered by SPO2_MARGIN_SAMPLES).
 *
 * @param in    input signal, length n
 * @param out   filtered signal, length n (may be the same buffer as in)
 * @param n     number of samples
 */
static void bandpass_filter(const float *in, float *out, int n)
{
    /* Delay-line state for each biquad section: w[s][0]/w[s][1] hold the
     * two feedback registers of the transposed Direct Form II structure. */
    float w[N_BIQUAD][2];
    /* Zero initial state: valid because the filter's startup transient
     * falls inside SPO2_MARGIN_SAMPLES and gets cropped out by the caller. */
    memset(w, 0, sizeof(w));

    /* One sample at a time, run through the full cascade of sections. */
    for (int i = 0; i < n; i++)
    {
        float x = in[i];
        /* Cascade the 4 biquad sections: each one's output feeds the next
         * one's input, together implementing the 4th-order response. */
        for (int s = 0; s < N_BIQUAD; s++)
        {
            /* Numerator coefficients of this section. */
            float b0 = SOS[s][0], b1 = SOS[s][1], b2 = SOS[s][2];
            /* Denominator coefficients of this section (a0 is always 1). */
            float a1 = SOS[s][4], a2 = SOS[s][5];
            /* Transposed Direct Form II: one multiply-add per state
             * register, avoiding the need to keep past input/output
             * samples explicitly. */
            float y = b0 * x + w[s][0];
            w[s][0] = b1 * x - a1 * y + w[s][1];
            w[s][1] = b2 * x - a2 * y;
            /* This section's output becomes the next section's input. */
            x = y;
        }
        out[i] = x;
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * 2. SINGLE-POINT DFT  (Goertzel over the cardiac band)
 *
 * Instead of a full FFT, we compute the DFT only at the frequencies of
 * interest (0.5-4.0 Hz in 0.1 Hz steps). This is O(N*K) with K≈35,
 * efficient on a microcontroller without an external FFT library.
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * Goertzel algorithm: spectral power at frequency f.
 * Includes a Hanning window to reduce spectral leakage.
 *
 * @param x     signal (with DC already removed)
 * @param n     length
 * @param f     target frequency (Hz)
 * @param fs    sampling frequency (Hz)
 * @return      power (squared units)
 */
static float goertzel_power(const float *x, int n, float f, float fs)
{
    /* Angular frequency of the target bin, and the Goertzel recurrence
     * coefficient derived from it. */
    float omega = 2.0f * (float)M_PI * f / fs;
    float coeff = 2.0f * cosf(omega);
    /* Recursion state: the two most recent intermediate values. */
    float s_prev = 0.0f, s_prev2 = 0.0f;

    for (int i = 0; i < n; i++)
    {
        /* Hanning window, applied sample-by-sample as it's fed into the
         * recursion (equivalent to windowing the whole signal first). */
        float w = 0.5f * (1.0f - cosf(2.0f * (float)M_PI * (float)i / (float)(n - 1)));
        float s = w * x[i] + coeff * s_prev - s_prev2;
        s_prev2 = s_prev;
        s_prev = s;
    }
    /* Closed-form Goertzel power from the final recursion state, divided
     * by n^2 to normalize away the signal length. */
    float power = (s_prev2 * s_prev2 + s_prev * s_prev - coeff * s_prev * s_prev2) / (float)(n * n);
    return power;
}

/**
 * Finds the dominant heart-rate frequency (Hz) within [FC_MIN, FC_MAX].
 * Resolution: 0.1 Hz.
 */
static float find_dominant_freq(const float *ac_signal, int n, float fs)
{
    float best_power = -1.0f;
    float best_freq = SPO2_FC_MIN;

    /* Scan the whole cardiac band one 0.1 Hz bin at a time, keeping the
     * frequency with the strongest response (i.e. the heart-rate estimate). */
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
 * RMS amplitude in the band [f_c - BW, f_c + BW].
 * Sums the power across the bins in the window and returns its sqrt.
 */
static float dft_amplitude(const float *signal_no_dc, int n, float fs,
                           float f_c, float bw)
{
    float power_sum = 0.0f;
    int count = 0;

    for (float f = f_c - bw; f <= f_c + bw + 0.01f; f += 0.1f)
    {
        /* Skip the near-zero-frequency bin: not meaningful once DC has
         * already been removed from the signal, and avoids a degenerate
         * bin when f_c - bw drops to (or below) 0 Hz. */
        if (f < 0.01f)
            continue;
        power_sum += goertzel_power(signal_no_dc, n, f, fs);
        count++;
    }
    /* Guard against an empty window (no valid bins): return 0 rather than
     * sqrt of an accumulator that was never touched. */
    return (count > 0) ? sqrtf(power_sum) : 0.0f;
}

/**
 * SNR (dB) of the cardiac peak against the band's noise floor, and
 * spectral purity: the fraction of the whole cardiac band's energy that's
 * concentrated in the peak window (1.0 = all the energy at the peak,
 * lower values = energy spread into broadband noise / motion).
 */
static void compute_snr_purity(const float *ac_signal, int n, float fs,
                               float f_c, float bw,
                               float *out_snr_db, float *out_purity)
{
    /* Strongest single bin inside the peak window (SNR numerator). */
    float peak_power = 0.0f;
    /* Sum of power across every bin inside the peak window (purity
     * numerator: how much energy actually sits at the detected heart rate). */
    float band_power = 0.0f;
    /* Sum of power across every bin in the full cardiac band (used both
     * as the SNR noise-floor estimate and as the purity denominator). */
    float total_power = 0.0f;
    int total_bins = 0;

    for (float f = SPO2_FC_MIN; f <= SPO2_FC_MAX + 0.01f; f += 0.1f)
    {
        float p = goertzel_power(ac_signal, n, f, fs);
        total_power += p;
        total_bins++;
        if (f >= f_c - bw && f <= f_c + bw)
        {
            band_power += p;
            if (p > peak_power)
                peak_power = p;
        }
    }

    /* Average power per bin across the whole band, used as the noise
     * floor: SNR compares the single strongest peak bin against this
     * average, not against the bins outside the peak window specifically. */
    float noise_mean = (total_bins > 0)
                           ? total_power / (float)total_bins
                           : 1e-12f;

    /* Guard against a degenerate (near-zero) noise floor before dividing. */
    *out_snr_db = (noise_mean < 1e-30f)
                      ? 0.0f
                      : 10.0f * log10f(peak_power / noise_mean);
    /* Guard against a silent signal (total_power ~ 0) before dividing. */
    *out_purity = (total_power > 1e-30f) ? (band_power / total_power) : 0.0f;
}

/* ═══════════════════════════════════════════════════════════════════════
 * 2b. SQI (SIGNAL QUALITY INDEX)
 * ═══════════════════════════════════════════════════════════════════════ */

/* Clamps x to [0, 1]; used to keep every SQI component well-formed before
 * it's blended and scaled to 8 bits. */
static float clamp01(float x)
{
    if (x < 0.0f)
        return 0.0f;
    if (x > 1.0f)
        return 1.0f;
    return x;
}

/**
 * Trapezoid membership function: 0 below low_min and above high_max,
 * linear ramp on the two transition segments, 1 on the flat plateau
 * [low_opt, high_opt]. Used to score a metric that's "bad" both when too
 * low and when too high (e.g. perfusion index).
 */
static float trapezoid(float x, float low_min, float low_opt,
                       float high_opt, float high_max)
{
    /* Outside the whole supported range: score is 0. */
    if (x <= low_min || x >= high_max)
        return 0.0f;
    /* Rising edge: between low_min and low_opt. */
    if (x < low_opt)
        return (x - low_min) / (low_opt - low_min);
    /* Falling edge: between high_opt and high_max. */
    if (x > high_opt)
        return (high_max - x) / (high_max - high_opt);
    /* Flat plateau: between low_opt and high_opt. */
    return 1.0f;
}

/**
 * Blends SNR, perfusion index and spectral purity into a 0-255 score.
 */
static uint8_t compute_sqi(float snr_db, float pi_pct, float purity)
{
    /* SNR component: linear ramp between SQI_SNR_MIN and SQI_SNR_MAX. */
    float f_snr = clamp01((snr_db - SQI_SNR_MIN) / (SQI_SNR_MAX - SQI_SNR_MIN));
    /* PI component: trapezoid, penalizing both weak and abnormally high PI. */
    float f_pi = trapezoid(pi_pct, SQI_PI_LOW_MIN, SQI_PI_LOW_OPT,
                           SQI_PI_HIGH_OPT, SQI_PI_HIGH_MAX);
    /* Spectral-purity component: already in [0,1] by construction, clamp01
     * here only guards against float rounding at the edges. */
    float f_purity = clamp01(purity);

    /* Weighted blend (weights defined in spo2_pipeline.h, sum to 1.0),
     * then clamp once more and scale from [0,1] to the 0-255 output range. */
    float sqi = SQI_W_SNR * f_snr + SQI_W_PI * f_pi + SQI_W_PURITY * f_purity;
    return (uint8_t)lroundf(clamp01(sqi) * 255.0f);
}

/* ═══════════════════════════════════════════════════════════════════════
 * 3. RATIO R COMPUTATION
 * ═══════════════════════════════════════════════════════════════════════ */

Spo2Status spo2_compute_R(const float *ir, const float *red, int n_samples,
                          float *out_R, float *out_snr, float *out_pi,
                          float *out_hr_bpm, uint8_t *out_sqi)
{
    /* Static buffers: no malloc, same embedded-friendly convention used
     * throughout this pipeline (deterministic, fixed memory footprint).
     * Sized to the maximum supported capture; the actually-used length is
     * n, derived below from n_samples (the real length of the raw signal). */
    static float ir_work[SPO2_MAX_N_SAMPLES];
    static float red_work[SPO2_MAX_N_SAMPLES];
    static float ir_ac[SPO2_MAX_N_SAMPLES];
    static float ir_nodC[SPO2_MAX_N_SAMPLES];
    static float red_nodC[SPO2_MAX_N_SAMPLES];

    const int m = SPO2_MARGIN_SAMPLES;
    /* Working window length: n_samples minus one margin on each side.
     * Replaces the old compile-time SPO2_N_WORK now that capture length
     * varies per raw signal instead of being fixed at 600. */
    const int n = n_samples - 2 * m;

    /* Reject captures too short to leave a usable window after cropping,
     * or longer than the statically-sized buffers can hold. */
    if (n_samples > SPO2_MAX_N_SAMPLES || n <= 1)
    {
        *out_R = 0.0f;
        *out_snr = 0.0f;
        *out_pi = 0.0f;
        *out_sqi = 0;
        return SPO2_ERR_N_SAMPLES;
    }

    /* Copy the working region, cropping SPO2_MARGIN_SAMPLES off each edge
     * of the raw capture (discards the bandpass filter's startup
     * transient and any edge artifacts from the acquisition window). */
    memcpy(ir_work, ir + m, (size_t)n * sizeof(float));
    memcpy(red_work, red + m, (size_t)n * sizeof(float));

    /* DC = arithmetic mean of the cropped window, per channel. */
    float dc_ir = 0.0f, dc_red = 0.0f;
    for (int i = 0; i < n; i++)
    {
        dc_ir += ir_work[i];
        dc_red += red_work[i];
    }
    dc_ir /= (float)n;
    dc_red /= (float)n;

    /* DC-removed signal for both channels, ready for spectral analysis. */
    for (int i = 0; i < n; i++)
    {
        ir_nodC[i] = ir_work[i] - dc_ir;
        red_nodC[i] = red_work[i] - dc_red;
    }

    /* Bandpass-filter the IR channel to isolate the cardiac band before
     * searching for the dominant frequency, so out-of-band noise/motion
     * doesn't distort the peak search. */
    bandpass_filter(ir_nodC, ir_ac, n);

    /* Dominant heart-rate frequency, from the filtered IR signal. */
    float f_c = find_dominant_freq(ir_ac, n, SPO2_FS);
    *out_hr_bpm = roundf(f_c * 60.0f); /* Hz -> beats per minute (integer) */

    /* AC amplitude via narrow-band DFT around f_c, for both channels
     * (dft_amplitude acts as its own narrowband filter, so this runs on
     * the unfiltered DC-removed signal directly). */
    float amp_ir = dft_amplitude(ir_nodC, n, SPO2_FS, f_c, SPO2_DFT_BW);
    float amp_red = dft_amplitude(red_nodC, n, SPO2_FS, f_c, SPO2_DFT_BW);

    if (amp_ir < 1e-6f)
    {
        /* No detectable pulsatile signal on IR: any ratio computed from
         * it would be meaningless, so zero out every output and bail. */
        *out_R = 0.0f;
        *out_snr = 0.0f;
        *out_pi = 0.0f;
        *out_sqi = 0;
        return SPO2_ERR_LOW_AC;
    }

    /* Ratio of ratios: the standard pulse-oximetry R value. */
    float R = (amp_red / dc_red) / (amp_ir / dc_ir);
    float snr, purity;
    /* SNR/purity computed from the filtered IR signal (ir_ac), the same
     * one used for the frequency search, so both metrics reflect how much
     * of that already-band-limited signal's power truly sits at f_c. */
    compute_snr_purity(ir_ac, n, SPO2_FS, f_c, SPO2_DFT_BW, &snr, &purity);
    /* IR perfusion index: AC/DC ratio of the infrared channel, as a percentage. */
    float pi = (amp_ir / dc_ir) * 100.0f;

    *out_R = R;
    *out_snr = snr;
    *out_pi = pi;
    *out_sqi = compute_sqi(snr, pi, purity);

    /* Quality gate: reject on low SNR or on R falling outside the
     * physiologically plausible range. Checked (and returned) in this
     * order, but every output above is already filled in either case. */
    if (snr < SPO2_SNR_THRESHOLD)
        return SPO2_ERR_LOW_SNR;
    if (R < SPO2_R_MIN || R > SPO2_R_MAX)
        return SPO2_ERR_R_RANGE;

    return SPO2_OK;
}

/* ═══════════════════════════════════════════════════════════════════════
 * 4. SpO2 PREDICTION
 * ═══════════════════════════════════════════════════════════════════════ */

float spo2_predict(float R)
{
    float spo2 = SPO2_CAL_0 + SPO2_CAL_1 * R + SPO2_CAL_2 * R * R;
    return fmaxf(SPO2_OUT_MIN, fminf(SPO2_OUT_MAX, spo2));
}

Spo2Status spo2_compute(const float *ir, const float *red, int n_samples,
                        float *out_spo2, float *out_hr_bpm, uint8_t *out_sqi)
{
    float R, snr, pi;
    Spo2Status status = spo2_compute_R(ir, red, n_samples, &R, &snr, &pi,
                                       out_hr_bpm, out_sqi);
    /* Only predict SpO2 once R has passed every quality check; on error,
     * *out_spo2 is deliberately left untouched (see header docstring). */
    if (status == SPO2_OK)
    {
        *out_spo2 = spo2_predict(R);
    }
    return status;
}
