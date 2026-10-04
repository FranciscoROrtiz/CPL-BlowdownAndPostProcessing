# filename: PP_ColdFlow2_PostProcessing_v8.py
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
file = r"C:\Users\franc\Downloads\TCF4.parquet"

# Name of directory to save figures of data in
out_dir = r"C:\Users\franc\Downloads\THF1"
os.makedirs(out_dir, exist_ok=True)

# Name of test for labeling
test_name = "Toad - HF1"

t_min = 0
t_max = 1200

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

# ============================================================
# PT PLOTS
# ============================================================
if Fuel_Run == True:
    start_idx, end_idx = MBV_start_end_idx(df['MBV_fuel_run.open'].values)
    plot_pt(df['time_sec'], 'PT_fuel_tank', ['SV_fuel_vent.open', 'MBV_fuel_run.open', 'SV_N2_fill.open'], 'PT_fuel_tank - Fuel_Run', [(start_idx/10)-1, (end_idx/10)+1])
else:
    plot_pt(df['time_sec'], 'PT_fuel_tank', ['SV_fuel_vent.open', 'MBV_fuel_run.open', 'SV_N2_fill.open'], 'PT_fuel_tank', [t_min, t_max])

if Ox_Run == True:
    start_idx, end_idx = MBV_start_end_idx(df['MBV_ox_run.open'].values)
    plot_pt(df['time_sec'], 'PT_ox_tank', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'PT_ox_tank - Ox_Run', [(start_idx/10)-1, (end_idx/10)+1])
    plot_pt(df['time_sec'], 'LC_ox_tank', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'LC_ox_tank - Ox_Run', [(start_idx/10)-1, (end_idx/10)+1])
else:
    plot_pt(df['time_sec'], 'PT_ox_tank', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'PT_ox_tank', [t_min, t_max])

if Fuel_Fill == True:
    start_idx, end_idx = MBV_start_end_idx(df['SV_N2_fill.open'].values)
    plot_pt(df['time_sec'], 'PT_n2_fill', ['SV_N2_fill.open', 'MBV_fuel_run.open', 'SV_fuel_vent.open'], 'PT_n2_fill - Fuel_Fill', [(start_idx/10)-1, (end_idx/10)+1])
else:
    plot_pt(df['time_sec'], 'PT_n2_fill', ['SV_N2_fill.open', 'MBV_fuel_run.open', 'SV_fuel_vent.open'], 'PT_n2_fill', [t_min, t_max])

if Ox_Fill == True:
    start_idx, end_idx = MBV_start_end_idx(df['SV_ox_fill.open'].values)
    plot_pt(df['time_sec'], 'PT_ox_fill', ['SV_ox_fill.open'], 'PT_ox_fill - Ox_Fill', [(start_idx/10)-1, (end_idx/10)+1])
else:
    plot_pt(df['time_sec'], 'PT_ox_fill', ['SV_ox_fill.open'], 'PT_ox_fill', [t_min, t_max])

plot_pt(df['time_sec'], 'PT_ox_fill', ['SV_ox_fill.open'], 'PT_ox_fill', [t_min, t_max])
plot_pt(df['time_sec'], 'PT_n2_fill', ['SV_N2_fill.open', 'MBV_fuel_run.open', 'SV_fuel_vent.open'], 'PT_n2_fill', [t_min, t_max])
plot_pt(df['time_sec'], 'PT_ox_tank', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'PT_ox_tank', [t_min, t_max])
plot_pt(df['time_sec'], 'PT_fuel_tank', ['SV_fuel_vent.open', 'MBV_fuel_run.open', 'SV_N2_fill.open'], 'PT_fuel_tank', [t_min, t_max])
plot_pt(df['time_sec'], 'PT_chamber', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'PT_chamber', [t_min, t_max])
plot_pt(df['time_sec'], 'LC_ox_tank', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'LC_ox_tank', [t_min, t_max])
plot_pt(df['time_sec'], 'TC_aft', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'TC_aft', [t_min, t_max])
plot_pt(df['time_sec'], 'TC_mid', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'TC_mid', [t_min, t_max])
plot_pt(df['time_sec'], 'TC_fwd', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'TC_fwd', [t_min, t_max])
plot_pt(df['time_sec'], 'pt_VT_throat', ['SV_fuel_vent.open', 'MBV_fuel_run.open', 'SV_N2_fill.open'], 'pt_VT_throat', [t_min, t_max])
plot_pt(df['time_sec'], 'pt_VT_upstream', ['SV_fuel_vent.open', 'MBV_fuel_run.open', 'SV_N2_fill.open'], 'pt_VT_upstream', [t_min, t_max])
plot_pt(df['time_sec'], 'LC_TOTAL', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'LC_TOTAL', [t_min, t_max])
# plt.show()
################################################################################

################################################################################
# BUTTER
df['Butter_PT_ox_tank'] = lowpass_filter(df['PT_ox_tank'], 1, 10)
df['Butter_PT_fuel_tank'] = lowpass_filter(df['PT_fuel_tank'], 1, 10)
df['Butter_PT_chamber'] = lowpass_filter(df['PT_chamber'], 2, 10)
df['Butter_PT_n2_fill'] = lowpass_filter(df['PT_n2_fill'], 3, 10)
df['Butter_LC_ox_tank'] = lowpass_filter(df['LC_ox_tank'], 0.5, 10)
df['Butter_LC_TOTAL'] = lowpass_filter(df['LC_TOTAL'], 0.5, 4)

if Butter == True:
    plot_pt(df['time_sec'], 'Butter_PT_ox_tank', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'Butter_PT_ox_tank')
    plot_pt(df['time_sec'], 'Butter_PT_fuel_tank', ['SV_fuel_vent.open', 'MBV_fuel_run.open', 'SV_N2_fill.open'], 'Butter_PT_fuel_tank')
    plot_pt(df['time_sec'], 'Butter_PT_chamber', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'Butter_PT_chamber')
    plot_pt(df['time_sec'], 'Butter_PT_n2_fill', ['SV_N2_fill.open', 'MBV_fuel_run.open', 'SV_fuel_vent.open'], 'Butter_PT_n2_fill')
    plot_pt(df['time_sec'], 'Butter_LC_ox_tank', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'Butter_LC_ox_tank')
    plot_pt(df['time_sec'], 'Butter_LC_TOTAL', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'Butter_LC_TOTAL')

# SMA
# df['sma'] = df['Butter_PT_ox_tank'].rolling(window=3).mean()
# plot_pt(df['time_sec'], 'sma', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'SMA3_PT_ox_tank')
# df['sma_5'] = df['Butter_PT_ox_tank'].rolling(window=5).mean()
# plot_pt(df['time_sec'], 'sma_5', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'SMA5_PT_ox_tank')


# ============================================================
# PSD PlOTS
# ============================================================
Fs = 10

if PSD == True:
    plt.figure()
    plt.psd(df['LC_ox_tank'], Fs=Fs)
    plt.title('PSD - LC_ox_tank')
    # plt.show()

    plt.figure()
    plt.psd(df['PT_n2_fill'], Fs=Fs)
    plt.title('PSD - PT_n2_fill')
    # plt.show()

    plt.figure()
    plt.psd(df['PT_chamber'], Fs=Fs)
    plt.title('PSD - PT_chamber')
    # plt.show()

    plt.figure()
    plt.psd(df['PT_ox_tank'], Fs=Fs)
    plt.title('PSD - PT_ox_tank')
    # plt.show()

    plt.figure()
    plt.psd(df['PT_fuel_tank'], Fs=Fs)
    plt.title('PSD - PT_fuel_tank')
    # plt.show()



# ============================================================
# LOAD CELLS
# ============================================================
# for lc in ['LC_chamber_1', 'LC_chamber_2', 'LC_chamber_3']:
#     if lc in df.columns:
#         plt.figure(figsize=(10,4))
#         plt.plot(df['time_sec'], df[lc])
#         plt.title(f"{test_name} - {lc}")
#         plt.grid()
#
#         apply_time_window()
#
#         save_fig(lc)





# ============================================================
# SHOW
# ============================================================
plt.show()
