# filename: Venturi_Calc.py
# created: 05.09.26
# last edited: 05.09.26
# authors: Francisco Ortiz, Diego Ortiz
# purpose: plot test data

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import CoolProp.CoolProp as CP
import os
import glob
from scipy.signal import savgol_filter

from scipy import signal

def lowpass_filter(data, cutoff, fs, order=15):
    # 'sos' output is recommended for better numerical stability
    sos = signal.butter(order, cutoff, btype='low', fs=fs, output='sos')
    # Use sosfiltfilt for zero-phase filtering (no time delay)
    filtered_data = signal.sosfiltfilt(sos, data)
    return filtered_data

def MBV_start_end_idx(MBV):
    MBV_idx = np.where(MBV == True)[0]
    start_idx = MBV_idx[0]
    end_idx = np.where((MBV == False) & (np.arange(len(MBV)) > start_idx))[0][0]
    return start_idx, end_idx
# ============================================================
# CONFIG
# ============================================================





"""
vvv FILE INPUTS vvv
"""
# Path to test data file
file = r"C:\Users\franc\Downloads\THF1.parquet"

# Name of directory to save figures of data in
out_dir = r"C:\Users\franc\Downloads\THF1"
os.makedirs(out_dir, exist_ok=True)

# Name of test for labeling
test_name = "Toad - HF1"

t_min = 520
t_max = 535

# If one is set to true, then graph will be generated according to the opening and closing of associated valve
Fuel_Run = True
Ox_Run = True
Fuel_Fill = False
Ox_Fill = False


# Under Maintenance # KEEP FALSE
PSD = False         # Set True if you want to display PSD of data
Butter = False      # Set True if you want to display data w/ Butterworth filter

"""
^^^ FILE INPUTS ^^^
"""





psi_to_pa = 6894.75729

# ============================================================
# CLEAN OUTPUT
# ============================================================
for f in glob.glob(os.path.join(out_dir, "*.png")):
    os.remove(f)

print("Deleted old plots")

# ============================================================
# LOAD DATA
# ============================================================
df = pd.read_parquet(file)







print("\nCOLUMNS IN DATAFRAME:")
for col in df.columns:
    print(col)

df['timestamp'] = pd.to_datetime(df['timestamp'])
df = df.sort_values('timestamp')
df['time_sec'] = (df['timestamp'] - df['timestamp'].iloc[0]).dt.total_seconds()
dt_sample = df['time_sec'].iloc[1]-df['time_sec'].iloc[0]
print(f"dt = {dt_sample} [s]")
print(f"f_s = {1/dt_sample} [Hz]")

# ============================================================
# SAVE FUNCTION
# ============================================================
def save_fig(name):
    path = os.path.join(out_dir, f"{test_name} - {name}.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')

# ============================================================
# TIME WINDOW HELPER
# ============================================================
def apply_time_window(t_min, t_max, ax=None):
    if ax is None:
        plt.xlim(t_min, t_max)
    else:
        ax.set_xlim(t_min, t_max)

# ============================================================
# PT PLOTS
# ============================================================
def plot_pt(time, pt_col, valve_cols, title, t_ends):
    t_min = t_ends[0]
    t_max = t_ends[1]
    fig, ax1 = plt.subplots(figsize=(10, 4))

    l1, = ax1.plot(time, df[pt_col], label=pt_col)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel(pt_col)
    ax1.grid(True)

    apply_time_window(t_min, t_max, ax1)

    ax2 = ax1.twinx()

    valve_lines = []
    valve_labels = []

    for v in valve_cols:
        if v in df.columns:
            line, = ax2.step(time, df[v], where='post', linestyle='--')
            valve_lines.append(line)
            valve_labels.append(v)

    ax2.set_ylabel("Valve State")

    ax1.legend([l1] + valve_lines, [pt_col] + valve_labels)

    plt.title(title)
    plt.tight_layout()
    # plt.xlim(210,230)
    save_fig(pt_col)


# idx_start = 5239
# idx_end = 5281
idx_start = 3765
idx_end = 3970
rho = 786
A1 = 0.0799229025 * 0.00064516
A2 = 0.02035830579 * 0.00064516
dP_VT = (df['pt_VT_upstream'].values[idx_start:idx_end] - df['pt_VT_throat'].values[idx_start:idx_end]) * 6894.76
sqrt_dP_VT = np.sqrt(dP_VT)
term = math.sqrt(2*rho) * sqrt_dP_VT
bottom = 1/(1-(A2/A1)**2)**0.5
mdot = (A2/bottom) * term
mdot_corr = (A2/bottom) * term * 0.8


plt.figure()
plt.plot(df['time_sec'].values[idx_start:idx_end], mdot)
plt.plot(df['time_sec'].values[idx_start:idx_end], mdot_corr)






plt.show()
