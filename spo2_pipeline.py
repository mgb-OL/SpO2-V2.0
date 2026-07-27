"""
Pipeline de cálculo de SpO2 a partir de señales PPG (IR + Red).

Método: extracción de AC vía FFT en banda estrecha alrededor de la
frecuencia cardíaca dominante. Más robusto que RMS broadband cuando
el perfusion index es bajo (<1%).

Parámetros validados con dataset de calibración (673 capturas, 70–99% SpO2):
  - R² = 0.72  (con umbral SNR ≥ 9 dB)
  - RMSE = 3.28%
  - Curva: SpO2 = 102.32 − 29.80 × R
"""

import re
import warnings

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, welch
from scipy.stats import linregress

warnings.filterwarnings("ignore")

# ──────────────────────────────────────────────
# Constantes
# ──────────────────────────────────────────────
FS = 50           # Hz — frecuencia de muestreo
MARGIN_S = 0.5    # segundos a descartar en los bordes tras filtrar
FC_MIN = 0.5      # Hz — límite inferior banda cardíaca
FC_MAX = 4.0      # Hz — límite superior banda cardíaca (~240 bpm)
FFT_BW = 0.1      # Hz — semiancho de ventana FFT alrededor de FC
SNR_THRESHOLD = 9 # dB — capturas por debajo se descartan
R_MIN = 0.1       # límite fisiológico inferior de R
R_MAX = 2.0       # límite fisiológico superior de R


# ──────────────────────────────────────────────
# Filtros
# ──────────────────────────────────────────────
def bandpass(signal: np.ndarray, fs: float = FS,
             low: float = FC_MIN, high: float = FC_MAX,
             order: int = 4) -> np.ndarray:
    """Filtro Butterworth paso banda (fase cero)."""
    b, a = butter(order, [low / (fs / 2), high / (fs / 2)], btype="band")
    return filtfilt(b, a, signal)


# ──────────────────────────────────────────────
# Extracción de AC por FFT
# ──────────────────────────────────────────────
def dominant_cardiac_freq(signal_ac: np.ndarray, fs: float = FS) -> float:
    """Frecuencia cardíaca dominante (Hz) por Welch PSD."""
    freqs, psd = welch(signal_ac, fs=fs, nperseg=min(256, len(signal_ac) // 2))
    mask = (freqs >= FC_MIN) & (freqs <= FC_MAX)
    return float(freqs[mask][np.argmax(psd[mask])])


def fft_amplitude(signal: np.ndarray, fs: float, f_cardiac: float,
                  bw: float = FFT_BW) -> float:
    """
    Amplitud RMS de la componente pulsátil en [f_cardiac-bw, f_cardiac+bw].
    Usa ventana Hanning para reducir spectral leakage.
    """
    n = len(signal)
    fft_vals = np.fft.rfft(signal * np.hanning(n))
    freqs = np.fft.rfftfreq(n, 1.0 / fs)
    mask = (freqs >= f_cardiac - bw) & (freqs <= f_cardiac + bw)
    return float(np.sqrt(np.sum(np.abs(fft_vals[mask]) ** 2)) / n * 2)


def signal_snr(signal_ac: np.ndarray, fs: float, f_cardiac: float,
               bw: float = FFT_BW) -> float:
    """SNR (dB) del pico cardíaco respecto al ruido de fondo en la banda."""
    freqs, psd = welch(signal_ac, fs=fs, nperseg=min(256, len(signal_ac) // 2))
    mask_peak = (freqs >= f_cardiac - bw) & (freqs <= f_cardiac + bw)
    mask_all  = (freqs >= FC_MIN) & (freqs <= FC_MAX)
    peak  = psd[mask_peak].max()
    noise = np.mean(psd[mask_all])
    return float(10 * np.log10(peak / noise)) if noise > 0 else 0.0


# ──────────────────────────────────────────────
# Cálculo del ratio R
# ──────────────────────────────────────────────
def compute_R(ir: np.ndarray, red: np.ndarray,
              fs: float = FS, margin_s: float = MARGIN_S
              ) -> tuple[float, float, float]:
    """
    Calcula el ratio R = (AC_red/DC_red) / (AC_ir/DC_ir).

    Returns
    -------
    R       : ratio de ratios
    snr_db  : SNR de la señal IR en la frecuencia cardíaca (dB)
    pi      : perfusion index IR (%)
    """
    margin = int(margin_s * fs)
    ir_t  = ir[margin:-margin]
    red_t = red[margin:-margin]

    dc_ir  = float(np.mean(ir_t))
    dc_red = float(np.mean(red_t))

    # Detectar FC sobre IR filtrado
    ac_ir_filtered = bandpass(ir_t, fs)
    f_c = dominant_cardiac_freq(ac_ir_filtered, fs)

    # Amplitud AC via FFT en banda estrecha
    amp_ir  = fft_amplitude(ir_t  - dc_ir,  fs, f_c)
    amp_red = fft_amplitude(red_t - dc_red, fs, f_c)

    if amp_ir < 1e-6:
        return np.nan, 0.0, 0.0

    R   = (amp_red / dc_red) / (amp_ir / dc_ir)
    snr = signal_snr(ac_ir_filtered, fs, f_c)
    pi  = amp_ir / dc_ir * 100

    return R, snr, pi


# ──────────────────────────────────────────────
# Filtro de calidad
# ──────────────────────────────────────────────
def is_valid(R: float, snr: float,
             snr_th: float = SNR_THRESHOLD,
             r_min: float = R_MIN, r_max: float = R_MAX) -> bool:
    """True si la captura supera los criterios de calidad."""
    return (
        not np.isnan(R)
        and r_min <= R <= r_max
        and snr >= snr_th
    )


# ──────────────────────────────────────────────
# Calibración (regresión lineal)
# ──────────────────────────────────────────────
def calibrate(R_values: np.ndarray,
              spo2_ref: np.ndarray) -> tuple[float, float, float]:
    """
    Ajusta SpO2 = a + b * R por regresión lineal.

    Returns
    -------
    a, b   : coeficientes de la curva
    r2     : coeficiente de determinación
    """
    slope, intercept, r_val, _, _ = linregress(R_values, spo2_ref)
    return float(intercept), float(slope), float(r_val ** 2)


def predict_spo2(R: float, a: float, b: float) -> float:
    """Devuelve SpO2 estimado (%) a partir de R y los coeficientes de calibración."""
    return a + b * R


# ──────────────────────────────────────────────
# Parseo del formato del dataset (rd.csv)
# ──────────────────────────────────────────────
def _extract_arrays(cell_str: str) -> tuple[np.ndarray | None,
                                             np.ndarray | None,
                                             float | None]:
    """Parsea una celda del CSV que contiene un dict con IR, Red y Reference_SpO2."""
    def get_array(key):
        pattern = rf"'{key}':\s*array\(\[([\d.,\s\n-]+)\]"
        m = re.search(pattern, cell_str)
        if not m:
            return None
        vals = [float(v.strip())
                for v in m.group(1).replace("\n", "").split(",")
                if v.strip()]
        return np.array(vals)

    ref_m = re.search(r"'Reference_SpO2':\s*([\d.]+)", cell_str)
    ref   = float(ref_m.group(1)) if ref_m else None
    return get_array("IR"), get_array("Red"), ref


def load_dataset(csv_path: str) -> pd.DataFrame:
    """
    Carga rd.csv y devuelve un DataFrame con columnas:
    R, SNR, PI, SpO2_ref
    Solo incluye capturas que superan el filtro de calidad.
    """
    df_raw = pd.read_csv(csv_path)
    records = []

    for i in range(len(df_raw)):
        for j in range(len(df_raw.columns)):
            cell = df_raw.iloc[i, j]
            if not isinstance(cell, str):
                continue
            ir, red, ref = _extract_arrays(cell)
            if ir is None or red is None or ref is None or len(ir) < 100:
                continue

            R, snr, pi = compute_R(ir, red)

            records.append({
                "R":        R,
                "SNR":      snr,
                "PI":       pi,
                "SpO2_ref": ref,
                "valid":    is_valid(R, snr),
            })

    return pd.DataFrame(records)


# ──────────────────────────────────────────────
# Main: calibración y evaluación
# ──────────────────────────────────────────────
if __name__ == "__main__":
    import sys

    csv_path = sys.argv[1] if len(sys.argv) > 1 else "rd.csv"
    print(f"Cargando {csv_path} ...")

    df = load_dataset(csv_path)
    print(f"Capturas totales procesadas : {len(df)}")
    print(f"Capturas válidas (SNR≥{SNR_THRESHOLD}dB, R∈[{R_MIN},{R_MAX}]): "
          f"{df.valid.sum()}")

    df_valid = df[df.valid].copy()

    # Calibración
    a, b, r2 = calibrate(df_valid.R.values, df_valid.SpO2_ref.values)
    print(f"\nCurva de calibración: SpO2 = {a:.2f} + {b:.2f} × R")
    print(f"R² = {r2:.4f}")

    # Evaluación
    df_valid["SpO2_pred"] = predict_spo2(df_valid.R.values, a, b)
    df_valid["error"]     = df_valid.SpO2_pred - df_valid.SpO2_ref

    print("\nRMSE por nivel de SpO2:")
    for v in sorted(df_valid.SpO2_ref.unique()):
        s = df_valid[df_valid.SpO2_ref == v]
        rmse_v = float(np.sqrt(np.mean(s.error ** 2)))
        print(f"  {v:.0f}%  n={len(s):3d}  RMSE={rmse_v:.2f}%")
        
    rmse = float(np.sqrt(np.mean(df_valid.error ** 2)))
    mae  = float(np.mean(np.abs(df_valid.error)))
    bias = float(np.mean(df_valid.error))

    print(f"\nMétricas globales:")
    print(f"  RMSE = {rmse:.3f}%")
    print(f"  MAE  = {mae:.3f}%")
    print(f"  Bias = {bias:.3f}%")

    
