import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter, find_peaks

# ==========================================
# 1. MATHEMATICAL DEFINITIONS
# ==========================================
def gaussian(x, amp, cen, wid):
    """1D Gaussian profile. wid is the standard deviation (sigma)."""
    return amp * np.exp(-(x - cen)**2 / (2 * wid**2))

def multi_gaussian(x, *params):
    """Sum of multiple Gaussian peaks."""
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        amp, cen, wid = params[i:i+3]
        y += gaussian(x, amp, cen, wid)
    return y

FWHM_FACTOR = 2.3548

# ==========================================
# 2. PROLINE-SPECIFIC BOUNDS & CONSTRAINTS
# ==========================================


from pathlib import Path
import re

def extract_peptide_sequence(file_path: str) -> str | None:
    """
    Extracts a peptide sequence from an absolute file path.
    Assumes standard uppercase amino acid letters separated by non-letters (e.g., underscores).
    """
    # 1. Extract just the filename without the extension (the "stem")
    filename = Path(file_path).stem
    
    # 2. Define the regex pattern for the 20 standard amino acids.
    # Negative lookarounds (?<!...) and (?!...) ensure the sequence isn't buried 
    # inside a standard English word (e.g., "DataACDFile").
    pattern = re.compile(r'(?<![A-Za-z])([ACDEFGHIKLMNPQRSTVWY]+)(?![A-Za-z])')
    
    # 3. Find all matches in the filename
    matches = pattern.findall(filename)
    
    if not matches:
        return None
        
    # 4. Return the longest match
    # (Filters out isolated single characters like "A" in "Experiment_A_ACDEF_data")
    return max(matches, key=len)

# --- Examples ---
#paths = [
#    "/usr/home/data/exp_1_ACDEFGHIK_results.csv",
#    "C:\\Users\\Lab\\run_A_MYSQLE_test.txt",
#    "/var/log/MALDI_WYF_001.raw"
#]

#for p in paths:
#    seq = extract_peptide_sequence(p)
#    print(f"Path: {p}\nSequence: {seq}\n")

import sys
exam_seq = extract_peptide_sequence(sys.argv[1])

peak_config = {
    '1_Tail':        [1560, 1600, 20, 40, 'exclude',       0.3], 
    '2_ProAmide(C=O)':    [1610, 1630, 15,  30, 'denom_only',       1.5], # J C P 2011, 135, 234507peak position and width
    '4_Core_Beta':     [1615, 1630, 15, 25, 'beta_fraction', 1.0], # Pro Sci 2004, 13, 3314/JPCL 2014, 5, 1984
    '5_PPII_Coil':     [1640, 1655, 20, 30, 'denom_only',    1.0], # Biopolymers 1986, 25, 469; Biochim Biophys Acta 2007, 1767, 1073 
    '6_Beta_Turn':     [1655, 1670, 15, 20, 'denom_only',    0.8], # Biopolymers 1986, 25, 469
    '7_PrimaryAmide(C=O)': [1670, 1685, 10, 15, 'exclude',       0.5],  # Biochim Biophys Acta 2007, 1767, 1073
    '8_High_Beta':     [1685, 1695, 10, 20, 'beta_fraction', 0.4], # Biochemical J 2009, 421, 415;  Molecules 2020, 25, 2498 
} if "P" in exam_seq else {
    '1_Anchor':        [1560, 1600, 20, 40, 'exclude',       0.3],
    '2_PrimaryAmide(NH2)':    [1610, 1625, 8,  15, 'exclude',       0.8], # Biochim Biophys Acta 2007, 1767, 1073; width --> Biopolymers, 1990, 30, 1243
    '4_Core_Beta':     [1610, 1630, 15, 25, 'beta_fraction', 1.0],
    '5_PPII_Coil':     [1640, 1655, 20, 30, 'denom_only',    1.0],
    '6_Beta_Turn':     [1655, 1670, 15, 20, 'denom_only',    0.8],
    '7_PrimaryAmide(C=O)': [1670, 1685, 10, 15, 'exclude',       0.5],
    '8_High_Beta':     [1685, 1695, 10, 20, 'beta_fraction', 0.4]
}

# ==========================================
# 3. DATA LOADING & PRE-PROCESSING
# ==========================================
def load_and_prep_data(filepath):
    """Loads CSV, extracts 1580-1720 cm^-1 window, and applies baseline correction."""
    try:
        df = pd.read_csv(filepath)
        df.columns = df.columns.str.upper()
        x = np.array(df['WAVENUMBER'].values, dtype=np.float64)
        y = df['ABSORBANCE'].values
    except Exception as e:
        # Fallback to generate synthetic data if file not found (for testing)
        print(f"Could not load {filepath}. Generating synthetic test data...")
        x = np.linspace(1560, 1720, 500)
        y = multi_gaussian(x, 0.05, 1585, 12, 0.20, 1620, 5, 0.25, 1624, 6, 
                           0.45, 1629, 8, 0.30, 1646, 10, 0.25, 1662, 7, 
                           0.20, 1675, 5, 0.15, 1690, 6)
        y += np.random.normal(0, 0.002, size=x.size)

    mask = (x >= 1580) & (x <= 1720)
    x_fit = x[mask]
    y_fit = y[mask]
 
#    y_fit_corr = y_fit   
    slope = (y_fit[-1] - y_fit[0]) / (x_fit[-1] - x_fit[0])
    baseline = y_fit[0] + slope * (x_fit - x_fit[0])
    y_fit_corr = y_fit - baseline
    y_fit_corr = np.clip(y_fit_corr, a_min=0, a_max=None)
    
    return x_fit, y_fit_corr

# ==========================================
# 4. SECOND DERIVATIVE PEAK FINDING
# ==========================================
def identify_peaks_via_deriv2(x, y, window=15, poly=3):
    """
    Calculates the 2nd derivative using Savitzky-Golay and identifies peak centers.
    Returns the derivative array, the peak dictionary, and indices for plotting.
    """
    d2y = savgol_filter(y, window_length=window, polyorder=poly, deriv=2)
    inverted_d2y = -d2y
    
    noise_threshold = np.max(inverted_d2y) * 0.05 
    peak_indices, _ = find_peaks(inverted_d2y, prominence=noise_threshold)
    
    detected_centers = x[peak_indices]
    detected_intensities = inverted_d2y[peak_indices]
    empirical_peaks = dict(zip(detected_centers, detected_intensities))
    
    return d2y, empirical_peaks, peak_indices

# ==========================================
# 5. DYNAMIC OPTIMIZATION 
# ==========================================
def perform_deconvolution(x, y, empirical_peaks):
    lower_bounds = []
    upper_bounds = []
    p0 = []
    max_amp = np.max(y)
    
    for key, bounds in peak_config.items():
        c_min, c_max, fwhm_min, fwhm_max, role, weight = bounds
        
        valid_empirical_peaks = {
            wv: intensity for wv, intensity in empirical_peaks.items() 
            if c_min <= wv <= c_max
        }
        
        if valid_empirical_peaks:
            best_center_guess = max(valid_empirical_peaks, key=valid_empirical_peaks.get)
            print(f"[{key}] using empirical 2nd deriv center: {best_center_guess:.1f} cm^-1")
        else:
            best_center_guess = (c_min + c_max) / 2.0
            print(f"[{key}] no 2nd deriv peak found. Defaulting to: {best_center_guess:.1f} cm^-1")

        lower_bounds.extend([0.0, c_min, fwhm_min / FWHM_FACTOR])
        upper_bounds.extend([max_amp * 1.5, c_max, fwhm_max / FWHM_FACTOR])
        p0.extend([(max_amp / 3.0) * weight, best_center_guess, ((fwhm_min + fwhm_max) / 2.0) / FWHM_FACTOR])
        
    popt, _ = curve_fit(
        multi_gaussian, x, y, 
        p0=p0, bounds=(lower_bounds, upper_bounds), method='trf', maxfev=15000
    )
    
    return popt

# ==========================================
# 6. INTEGRATION & PLOTTING 
# ==========================================
def calculate_beta_fraction(popt):
    areas = {}
    total_backbone_area = 0.0
    beta_area = 0.0
    peak_names = list(peak_config.keys())
    
    for i in range(len(peak_names)):
        amp, cen, wid = popt[i*3 : i*3+3]
        area = amp * wid * np.sqrt(2 * np.pi)
        
        name = peak_names[i]
        areas[name] = area
        role = peak_config[name][4]
        
        if role in ['denom_only', 'beta_fraction']:
            total_backbone_area += area
        if role == 'beta_fraction':
            beta_area += area
                
    beta_fraction = (beta_area / total_backbone_area) * 100 if total_backbone_area > 0 else 0
    return areas, beta_fraction, total_backbone_area

def calculate_beta_orientation(areas, k=0.15):
    core_beta = areas.get('4_Core_Beta', areas.get('3_Core_Beta', 0.0))
    high_beta = areas.get('8_High_Beta', areas.get('7_High_Beta', 0.0))
    
    if core_beta == 0 and high_beta == 0: return 0.0, 0.0
        
    anti_low = high_beta / k
    total_anti = high_beta + anti_low
    total_para = core_beta - anti_low
    
    if total_para < 0:
        total_para = 0.0
        total_anti = core_beta + high_beta 
        
    total_assembled = total_anti + total_para
    pct_anti = (total_anti / total_assembled) * 100 if total_assembled > 0 else 0
    pct_para = (total_para / total_assembled) * 100 if total_assembled > 0 else 0
    
    return pct_anti, pct_para

def plot_results(x, y, popt, beta_fraction, pct_para, filepath=None):
    fig, ax = plt.subplots(figsize=(8, 6))

    y_fit = multi_gaussian(x, *popt)
    ax.plot(x, y, 'ko', markersize=9, label='Experimental Data', alpha=0.5)
    ax.plot(x, y_fit, 'r-', lw=4, label='Total Fit Envelope')

    peak_names = list(peak_config.keys())
    colors = plt.cm.viridis(np.linspace(0, 1, len(peak_names)))
    colors = ['sandybrown','gray','pink','skyblue','lightgreen','gray','pink'] 

    out_d={}
    out_d['wavenumer']= x
    out_d['ExptAbs'] = y
    out_d['fitAbs'] = y_fit
    for i in range(len(peak_names)):
        amp, cen, wid = popt[i*3 : i*3+3]
        y_peak = gaussian(x, amp, cen, wid)
        name = peak_names[i]
        out_d[name] = y_peak
    df = pd.DataFrame(out_d)
    df.to_csv(filepath[:-4]+"_FTIR_SourceData.csv", index=False)


    for i in range(len(peak_names)):
        amp, cen, wid = popt[i*3 : i*3+3]
        y_peak = gaussian(x, amp, cen, wid)
        name = peak_names[i]
        role = peak_config[name][4]
        
        if role == 'beta_fraction':
            ls, alpha, fill_alpha = '-', 0.9, 0.6
        elif role == 'denom_only':
            ls, alpha, fill_alpha = '--', 0.8, 0.4
        else:
            ls, alpha, fill_alpha = ':', 0.8, 0.3
            
        ax.fill_between(x, y_peak, alpha=fill_alpha, color=colors[i])
        ax.plot(x, y_peak, color=colors[i], lw=3, linestyle=ls, alpha=alpha, label=name[2:])
                 
    ax.set_xlim(1580, 1720)
    ax.set_xlabel('Wavenumber (cm$^{-1}$)', fontsize=32)
    ax.set_ylabel('Normalized Abs', fontsize=32)
    #ax.legend(loc='upper right', fontsize=18)

    for spine in ax.spines.values():
        spine.set_linewidth(4.0)
    ax.tick_params(axis='both', which='major', width=4.0, length=6, labelsize=32)
    
    plt.tight_layout()
    
    if filepath is None or filepath.endswith('test_data.csv'):
      plt.show()
    else:
      plt.savefig(filepath, dpi=100, bbox_inches='tight')

# ==========================================
# MAIN EXECUTION
# ==========================================
if __name__ == "__main__":
    fnm = sys.argv[1] if len(sys.argv) > 1 else "test_data.csv"
    
    x_data, y_data = load_and_prep_data(fnm)
    
    print("--- 2ND DERIVATIVE ANALYSIS ---")
    d2y_data, empirical_peaks, peak_indices = identify_peaks_via_deriv2(x_data, y_data, window=20, poly=3)
    
    print("\n--- OPTIMIZATION ---")
    optimized_params = perform_deconvolution(x_data, y_data, empirical_peaks)
    
    areas, beta_pct, tot_area = calculate_beta_fraction(optimized_params)
    pct_anti, pct_para = calculate_beta_orientation(areas, k=0.15)
    
    print("\n--- FIT RESULTS ---")
    peak_names = list(peak_config.keys())
    for i in range(len(peak_names)):
        name = peak_names[i]
        cen = optimized_params[i*3 + 1]
        fwhm = optimized_params[i*3 + 2] * FWHM_FACTOR
        role = peak_config[name][4].upper()
        print(f"{name:<16}: Center = {cen:6.1f} | FWHM = {fwhm:5.1f} | Area = {areas[name]:.3f} [{role}]")
    
    print("-" * 55)
    print(f"Corrected Backbone Beta-Sheet Content {fnm} : {beta_pct:.1f} % para {pct_para:.1f} %")
    print(f" -> Anti-parallel Sub-population  : {pct_anti:.1f} %")
    print(f" -> Parallel Sub-population       : {pct_para:.1f} %")
    print("-" * 55)
    
    plot_filepath = None if fnm == "test_data.csv" else fnm[:-4] + '_fit.png'
    plot_results(x_data, y_data, optimized_params, beta_pct, pct_para, filepath=plot_filepath)
