# filename: PP_ColdFlow2_PostProcessing_v8.py
# created: 05.09.26
# last edited: 05.09.26
# authors: Francisco Ortiz, Diego Ortiz
# purpose: plot test data

import json
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
file = r"C:\Users\franc\Downloads\CPL - Post Processing Scripts\All Coldflow and Hotfire Data\BG_HF_2.parquet"

# Name of directory to save figures of data in
out_dir = r"C:\Users\franc\Downloads\CPL - Post Processing Scripts\All Coldflow and Hotfire Data\BG_CF_1"
os.makedirs(out_dir, exist_ok=True)

# Name of test for labeling
test_name = f"BG - HF2"

t_min = 0
t_max = 5000

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


plot_pt(df['time_sec'], 'PT_ox_fill', ['SV_ox_fill.open'], 'PT_ox_fill', [t_min, t_max])
plot_pt(df['time_sec'], 'PT_n2_fill', ['SV_N2_fill.open', 'MBV_fuel_run.open', 'SV_fuel_vent.open'], 'PT_n2_fill', [t_min, t_max])
plot_pt(df['time_sec'], 'PT_ox_tank', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'PT_ox_tank', [t_min, t_max])
plot_pt(df['time_sec'], 'PT_fuel_tank', ['SV_fuel_vent.open', 'MBV_fuel_run.open', 'SV_N2_fill.open'], 'PT_fuel_tank', [t_min, t_max])
plot_pt(df['time_sec'], 'PT_chamber', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'PT_chamber', [t_min, t_max])
plot_pt(df['time_sec'], 'LC_ox_tank', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'LC_ox_tank', [t_min, t_max])
# plot_pt(df['time_sec'], 'TC_aft', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'TC_aft', [t_min, t_max])
# plot_pt(df['time_sec'], 'TC_mid', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'TC_mid', [t_min, t_max])
# plot_pt(df['time_sec'], 'TC_fwd', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'TC_fwd', [t_min, t_max])
# plot_pt(df['time_sec'], 'pt_VT_throat', ['SV_fuel_vent.open', 'MBV_fuel_run.open', 'SV_N2_fill.open'], 'pt_VT_throat', [t_min, t_max])
# plot_pt(df['time_sec'], 'pt_VT_upstream', ['SV_fuel_vent.open', 'MBV_fuel_run.open', 'SV_N2_fill.open'], 'pt_VT_upstream', [t_min, t_max])
plot_pt(df['time_sec'], 'LC_TOTAL', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'LC_TOTAL', [t_min, t_max])



# IMPORT SIMULATION DATA #############################################
# _TD_HF_1
# _BG_HF_1
json_file = f"BG_HF_2"
with open(f"t_{json_file}.json", 'r') as file:
    t_sim = json.load(file)
with open(f"Pc_{json_file}.json", 'r') as file:
    Pc_sim = json.load(file)
with open(f"Thrust_{json_file}.json", 'r') as file:
    Thrust_sim = json.load(file)
with open(f"Pox_{json_file}.json", 'r') as file:
    Pox_sim = json.load(file)
with open(f"Pfuel_{json_file}.json", 'r') as file:
    Pfuel_sim = json.load(file)
with open(f"mox_{json_file}.json", 'r') as file:
    mox_sim = json.load(file)
with open(f"Pc_{json_file}_adj.json", 'r') as file:
    Pc_sim_adj = json.load(file)
with open(f"Thrust_{json_file}_adj.json", 'r') as file:
    Thrust_sim_adj = json.load(file)

# # BG_HF_1
# t_shift = 755
# m_shift = -0.1777

# BG_HF_2
t_shift = 240
t_shift_2 = 1
m_shift = 0.433

# # BG_HF_3
# m_shift = 0.479
# t_shift = 216
# t_shift_2 = 1

# # TD_HF_1
# m_shift = 2.577
# t_shift = 523
# t_shift_2 = 1

time_sec = np.array(df['time_sec'])


sim_len = len(t_sim)

time_sec = time_sec - t_shift
idx_0 = np.where(time_sec >= t_sim[0])[0][0]
idx_f = np.where(time_sec >= t_sim[-1])[0][0]
print(f"{idx_0} & {idx_f}")
t_sim = np.array(t_sim)
# t_sim = np.array(t_sim) + t_shift
Pc_sim = np.array(Pc_sim)*(1/6894.76) - 14.7
Thrust_sim = np.array(Thrust_sim) * 0.224809
Pox_sim = np.array(Pox_sim)*(1/6894.76) - 14.7
Pfuel_sim = np.array(Pfuel_sim)*(1/6894.76) - 14.7
mox_sim = np.array(mox_sim)*2.20462 + m_shift

Pc_sim_adj = np.array(Pc_sim_adj)*(1/6894.76) - 14.7*np.ones(len(Pc_sim_adj))
Thrust_sim_adj = np.array(Thrust_sim_adj) * 0.224809


######################################################################



plt.figure()
plt.scatter(time_sec, df['PT_ox_tank'], label='data', c='k', marker='.')
plt.plot(t_sim, Pox_sim, label='sim', c='tab:red')
plt.title(f'Oxidizer Tank Pressure - {test_name}')
plt.xlabel('time [s]')
plt.ylabel('Pressure [psi]')
plt.legend()

plt.figure()
plt.scatter(time_sec, df['PT_fuel_tank'], label='data', c='k', marker='.')
plt.plot(t_sim, Pfuel_sim, label='sim', c='tab:red')
plt.title('Fuel Tank Pressure')
plt.xlabel('time [s]')
plt.ylabel('Pressure [psi]')
plt.legend()

# t_sim = t_sim + t_shift_2
time_sec = time_sec - t_shift_2

plt.figure()
plt.scatter(time_sec, df['PT_chamber'], label='data', c='k', marker='.')
plt.plot(t_sim, Pc_sim, label='sim', c='tab:red')
# plt.plot(t_sim, Pc_sim_adj, label='sim (n_c*=0.85)', c='tab:blue')
plt.title(f'Chamber Pressure - {test_name}')
plt.xlabel('time [s]')
plt.ylabel('Pressure [psi]')
plt.legend()

plt.figure()
plt.scatter(time_sec, df['LC_TOTAL'], label='data', c='k', marker='.')
plt.plot(t_sim, Thrust_sim, label='sim', c='tab:red')
# plt.plot(t_sim, Thrust_sim_adj, label='sim (n_c*=0.85)', c='tab:blue')
plt.title(f'Thrust - {test_name}')
plt.xlabel('time [s]')
plt.ylabel('Thrust [lb]')
plt.legend()

plt.figure()
plt.scatter(time_sec, df['LC_ox_tank'], label='data', c='k', marker='.')
plt.plot(t_sim, mox_sim, label='sim', c='tab:red')
plt.title(f'Oxidizer Mass - {test_name}')
plt.xlabel('time [s]')
plt.ylabel('Weight [lb]')
plt.legend()

print(df['time_sec'])


# ============================================================
# ERROR ANALYSIS
# ============================================================
T_sim_m = Thrust_sim[::10]
t_sim_m = t_sim[::10]
T_sim_um = T_sim_m[::5]
t_sim_um = t_sim_m[::5]
time_sec_um = time_sec[idx_0:idx_f+1:5]
T_um = np.array(df['LC_TOTAL'])[idx_0:idx_f+1:5]

T_sim_m_adj = Thrust_sim_adj[::10]
T_sim_um_adj = T_sim_m_adj[::5]

err_sim = (T_um-T_sim_um)/T_um
err_sim_adj = (T_um-T_sim_um_adj)/T_um

area_err = np.trapz(err_sim[5:], x=time_sec_um[5:])
err_avg = (1/(time_sec_um[-1]-time_sec_um[5]))*area_err
area_err_adj = np.trapz(err_sim_adj[5:], x=time_sec_um[5:])
err_avg_adj = (1/(time_sec_um[-1]-time_sec_um[5]))*area_err_adj

print(f"err={err_avg}")
print(f"err={err_avg_adj}")

plt.figure()
plt.scatter(time_sec_um, T_um)
plt.scatter(t_sim_um, T_sim_um)

plt.figure()
plt.scatter(time_sec_um, err_sim)
plt.scatter(time_sec_um, err_sim_adj, label='n_c* = 0.85')
plt.legend()

Pc_sim_m = Pc_sim[::10]
Pc_um = np.array(df['PT_chamber'])[idx_0:idx_f+1:5]


# ============================================================
# SHOW
# ============================================================
plt.show()
