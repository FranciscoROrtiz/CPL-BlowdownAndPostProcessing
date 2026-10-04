# filename: PP - Cold Flow 2 - Post Processing (1).py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import CoolProp.CoolProp as CP
import os
import glob
from scipy.signal import savgol_filter

# ============================================================
# CONFIG
# ============================================================
file = r"C:\Users\ddort\Downloads\Prop - Revel Parquets\Balrog-CF3.parquet"

out_dir = r"C:\Users\ddort\Downloads\Prop - Revel Parquets\Balrog\Cold Flow 1"
os.makedirs(out_dir, exist_ok=True)

test_name = "CF1"
psi_to_pa = 6894.75729

# GLOBAL TIME WINDOW
t_min = 1750
t_max = 2600

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

# ============================================================
# SAVE FUNCTION
# ============================================================
def save_fig(name):
    path = os.path.join(out_dir, f"Balrog - {test_name} - {name}.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')

# ============================================================
# TIME WINDOW HELPER
# ============================================================
def apply_time_window(ax=None):
    if ax is None:
        plt.xlim(t_min, t_max)
    else:
        ax.set_xlim(t_min, t_max)

# ============================================================
# PT PLOTS
# ============================================================
def plot_pt(time, pt_col, valve_cols, title):

    fig, ax1 = plt.subplots(figsize=(10, 4))

    l1, = ax1.plot(time, df[pt_col], label=pt_col)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel(pt_col)
    ax1.grid(True)

    apply_time_window(ax1)

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

    save_fig(pt_col)

# ============================================================
# PT PLOTS
# ============================================================
plot_pt(df['time_sec'], 'PT_ox_fill', ['SV_ox_fill.open'], 'PT_ox_fill')
plot_pt(df['time_sec'], 'PT_n2_fill', ['SV_N2_fill.open'], 'PT_n2_fill')
plot_pt(df['time_sec'], 'PT_chamber', ['MBV_ox_run.open', 'MBV_fuel_run.open'], 'PT_chamber')
plot_pt(df['time_sec'], 'PT_ox_tank', ['SV_ox_vent.open', 'MBV_ox_run.open', 'SV_ox_fill.open'], 'PT_ox_tank')
plot_pt(df['time_sec'], 'PT_fuel_tank', ['SV_fuel_vent.open', 'MBV_fuel_run.open', 'SV_N2_fill.open'], 'PT_fuel_tank')

# ============================================================
# LOAD CELLS
# ============================================================
for lc in ['LC_chamber_1', 'LC_chamber_2', 'LC_chamber_3']:
    if lc in df.columns:
        plt.figure(figsize=(10,4))
        plt.plot(df['time_sec'], df[lc])
        plt.title(f"{test_name} - {lc}")
        plt.grid()

        apply_time_window()

        save_fig(lc)

# ============================================================
# FUEL MASS FLOW (SPI)
# ============================================================
p1_fuel = df['PT_fuel_tank'].values * psi_to_pa
p2_fuel = df['PT_chamber'].values * psi_to_pa

Cd = 0.6
d_orifice = 0.001016
A_fuel = 3 * (math.pi / 4) * d_orifice**2

rho_fuel = CP.PropsSI('D', 'T', 298, 'P', np.mean(p1_fuel), "Water")

df['mdot_fuel_SPI'] = [
    0 if p1 <= p2 else Cd * A_fuel * math.sqrt(2 * rho_fuel * (p1 - p2))
    for p1, p2 in zip(p1_fuel, p2_fuel)
]

plt.figure()
plt.plot(df['time_sec'], df['mdot_fuel_SPI'])
plt.title(f"{test_name} - mdot_fuel")
plt.grid()

apply_time_window()

save_fig("mdot_fuel")

# ============================================================
# OX SPI
# ============================================================
fluid_ox = "REFPROP::CO2"
a0 = 0.00000281070638655

p1_ox = np.clip(df['PT_ox_tank'].values * psi_to_pa, 0, None)
p2_ox = np.clip(df['PT_chamber'].values * psi_to_pa, 0, None)

mdot_ox_SPI = []

for P1, P2 in zip(p1_ox, p2_ox):
    if P1 <= P2 or P1 < 1e5:
        mdot_ox_SPI.append(0)
        continue
    try:
        s = CP.PropsSI("S", "P", P1, "Q", 0, fluid_ox)
        rho1 = CP.PropsSI("D", "P", P1, "S", s, fluid_ox)
        mdot_ox_SPI.append(Cd * a0 * math.sqrt(2 * rho1 * (P1 - P2)))
    except:
        mdot_ox_SPI.append(0)

df['mdot_ox_SPI'] = mdot_ox_SPI

# ============================================================
# OX DYER
# ============================================================
mdot_ox_DYER = []

for P1, P2 in zip(p1_ox, p2_ox):
    if P1 <= P2 or P1 < 1e5:
        mdot_ox_DYER.append(0)
        continue
    try:
        T_sat = CP.PropsSI("T", "P", P1, "Q", 0, fluid_ox)
        s = CP.PropsSI("S", "P", P1, "Q", 0, fluid_ox)

        rho1 = CP.PropsSI("D", "P", P1, "S", s, fluid_ox)
        rho2 = CP.PropsSI("D", "P", P2, "S", s, fluid_ox)

        h1 = CP.PropsSI("H", "P", P1, "S", s, fluid_ox)
        h2 = CP.PropsSI("H", "P", P2, "S", s, fluid_ox)

        Pv = CP.PropsSI("P", "T", T_sat, "Q", 0, fluid_ox)

        deltaH = max(h1 - h2, 0.0)

        m_spi = Cd * a0 * math.sqrt(2 * rho1 * (P1 - P2))
        m_hem = Cd * a0 * rho2 * math.sqrt(2 * deltaH)

        k = math.sqrt((P1 - P2) / max(Pv - P2, 1e-6))

        mdot_ox_DYER.append((k/(1+k))*m_spi + (1/(1+k))*m_hem)

    except:
        mdot_ox_DYER.append(0)

df['mdot_ox_DYER'] = mdot_ox_DYER

# ============================================================
# OX COMPARISON
# ============================================================
plt.figure(figsize=(10, 4))

plt.plot(df['time_sec'], df['mdot_ox_SPI'],
         label='MDOT_OX_SPI', linewidth=2)

plt.plot(df['time_sec'], df['mdot_ox_DYER'],
         label='MDOT_OX_DYER', linestyle='--', linewidth=2)

plt.title(f"{test_name} - mdot_ox")
plt.xlabel("Time (s)")
plt.ylabel("Mass Flow Rate (kg/s)")
plt.grid(True)
plt.legend()

apply_time_window()

save_fig("mdot_ox")

# ============================================================
# OF RATIOS
# ============================================================
eps = 1e-9

df['OF_SPI'] = df['mdot_ox_SPI'] / (df['mdot_fuel_SPI'] + eps)
df['OF_DYER'] = df['mdot_ox_DYER'] / (df['mdot_fuel_SPI'] + eps)

# ============================================================
# OF PLOT
# ============================================================
plt.figure(figsize=(10, 4))

plt.plot(df['time_sec'], df['OF_SPI'],
         label='OF_SPI', linewidth=2)

plt.plot(df['time_sec'], df['OF_DYER'],
         label='OF_DYER', linestyle='--', linewidth=2)

plt.title(f"{test_name} - OF")
plt.xlabel("Time (s)")
plt.ylabel("OF")
plt.grid(True)
plt.legend()

apply_time_window()

save_fig("OF")

# ============================================================
# OX LOAD CELL → MASS FLOW RATE
# ============================================================

from scipy.signal import savgol_filter

# ------------------------------------------------------------
# 1. RAW LOAD CELL (lb → kg)
# ------------------------------------------------------------
m_ox_tank = df['LC_ox_tank'].values * 0.45359237
t = df['time_sec'].values

# ------------------------------------------------------------
# 2. SMOOTH SIGNAL (reduces derivative noise)
# ------------------------------------------------------------
m_ox_smooth = savgol_filter(
    m_ox_tank,
    window_length=31,   # must be odd; adjust if needed
    polyorder=3
)

# ------------------------------------------------------------
# 3. DIFFERENTIATE → MASS FLOW RATE (kg/s)
# ------------------------------------------------------------
dm_dt = np.gradient(m_ox_smooth, t)

# Positive outflow convention
df['MDOT_OX_TANK'] = -dm_dt


# ============================================================
# OX TANK MDOT (LOAD CELL) + VALVES OVERLAY
# ============================================================

fig, ax1 = plt.subplots(figsize=(10, 4))

# Main signal
l1, = ax1.plot(
    df['time_sec'],
    df['MDOT_OX_TANK'],
    label='MDOT_OX_TANK',
    linewidth=2
)

ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Mass Flow Rate (kg/s)")
ax1.grid(True)

apply_time_window(ax1)

# ------------------------------------------------------------
# Valve overlay (same style as PT_ox_tank)
# ------------------------------------------------------------
ax2 = ax1.twinx()

valve_cols = [
    'SV_ox_vent.open',
    'MBV_ox_run.open',
    'SV_ox_fill.open'
]

valve_lines = []
valve_labels = []

for v in valve_cols:
    if v in df.columns:
        line, = ax2.step(
            df['time_sec'],
            df[v],
            where='post',
            linestyle='--'
        )
        valve_lines.append(line)
        valve_labels.append(v)

ax2.set_ylabel("Valve State")

ax1.legend([l1] + valve_lines, ['MDOT_OX_TANK'] + valve_labels)

plt.title(f"{test_name} - MDOT_OX_TANK (Load Cell)")
plt.tight_layout()

save_fig("MDOT_OX_TANK")

# ============================================================
# SHOW
# ============================================================
plt.show()
