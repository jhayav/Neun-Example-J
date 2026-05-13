#!/usr/bin/env python3
"""
Estudio de Invariantes Dinámicos - CPG + CGC (Lymnaea stagnalis)
================================================================
Basado en: Garrido-Peña, Elices & Varona, Neurocomputing 461 (2021)

Ejecuta 4 condiciones de estimulación (N1M, N3t, SO, CGC) y genera
para cada una: box-plots de intervalos (Fig.5/7/9) y scatter plots
de correlación intervalo vs periodo con R (Fig.6/8/10).

Dependencias: numpy, pandas, matplotlib, seaborn, scipy
"""

import subprocess
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from scipy.signal import find_peaks
from pathlib import Path
import sys
import os


PROJECT_DIR = Path(__file__).resolve().parent
BUILD_DIR = PROJECT_DIR / "build"
OUTPUT_DIR = PROJECT_DIR / "resultados_invariantes"
OUTPUT_DIR.mkdir(exist_ok=True)

TARGETS = ['N1M', 'N3t', 'SO', 'CGC']

COL_T, COL_N1M, COL_N2V, COL_N3T, COL_SO, COL_CGC, COL_RAMP = range(7)

SPIKE_THRESHOLDS = {'N1M': -20.0, 'N2v': -25.0, 'N3t': -20.0}
MAX_ISI          = {'N1M': 300.0, 'N2v': 400.0, 'N3t': 600.0}
MIN_SPIKES       = {'N1M': 5,     'N2v': 1,     'N3t': 5}

TRANSIENT_TIME = 5000.0

RAMP_CONFIG = {
    'N1M': {
        'I_SO': -8.5,  'I_N1M': 'RAMPA', 'I_N2v': -2.0,  'I_N3t': 0.0,  'I_CGC': 0.03,
        'c_min': 0.0,  'c_max': -10.5,   'n_steps': 21,
        'nota': 'Paper'
    },
    'N3t': {
        'I_SO': -9.0,  'I_N1M': -10.0,   'I_N2v': -1.0,  'I_N3t': 'RAMPA', 'I_CGC': 0.03,
        'c_min': 0.0,  'c_max': -5.0,    'n_steps': 20,
        'nota': 'Paper'
    },
    'SO': {
        'I_SO': 'RAMPA', 'I_N1M': -10.0, 'I_N2v': -1.0,  'I_N3t': -4.0,  'I_CGC': 0.03,
        'c_min': -8.2, 'c_max': -13.0,   'n_steps': 20,
        'nota': 'Paper'
    },
    'CGC': {
        'I_SO': -8.5,  'I_N1M': -6.0,    'I_N2v': -2.0,  'I_N3t': 0.0,   'I_CGC': 'RAMPA',
        'c_min': 0.0,  'c_max': 0.15,    'n_steps': 20,
        'nota': 'Yo'
    },
}
RAMP_INTERVAL = 4600.0

sns.set_theme(style='whitegrid', context='paper', font_scale=1.0)
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 9,
    'axes.titlesize': 9,
    'axes.labelsize': 9,
    'figure.dpi': 150,
})

LABELS_PRETTY = {
    'period':     'N1M-Period',
    'N1M_BD':     'N1M-BD',   'N2v_BD':     'N2v-BD',   'N3t_BD':     'N3t-BD',
    'N1_N2_int':  'N1N2-Interval',    'N1_N3_int':  'N1N3-Interval',    'N2_N3_int':  'N2N3-Interval',
    'N2_N1_int':  'N2N1-Interval',  'N3_N1_int':  'N3N1-Interval',  'N3_N2_int':  'N3N2-Interval',
    'N1_N2_del':  'N1N2-Delay',  'N1_N3_del':  'N1N3-Delay',  'N2_N3_del':  'N2N3-Delay',
    'N2_N1_del':  'N2N1-Delay', 'N3_N1_del':  'N3N1-Delay', 'N3_N2_del':  'N3N2-Delay',
}

INTERVAL_KEYS = [
    'N1M_BD', 'N2v_BD', 'N3t_BD',
    'N1_N2_int', 'N1_N3_int', 'N2_N3_int',
    'N2_N1_int', 'N3_N1_int', 'N3_N2_int',
    'N1_N2_del', 'N1_N3_del', 'N2_N3_del',
    'N2_N1_del', 'N3_N1_del', 'N3_N2_del',
]

BOXPLOT_KEYS = [
    'period',
    'N1M_BD', 'N2v_BD', 'N3t_BD',
    'N1_N2_int', 'N1_N2_del',
    'N2_N1_int', 'N2_N1_del',
    'N1_N3_int', 'N1_N3_del',
    'N3_N1_int', 'N3_N1_del',
    'N2_N3_int', 'N2_N3_del',
    'N3_N2_int', 'N3_N2_del'
]


def compilar():
    """Compila el proyecto (cmake + make)."""
    print("[INFO] Compilando...")
    BUILD_DIR.mkdir(exist_ok=True)
    r = subprocess.run(["cmake", ".."], cwd=str(BUILD_DIR),
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("[ERROR] cmake falló:\n", r.stderr)
        sys.exit(1)
    r = subprocess.run(["make", "-j4", "cgc_ramp"], cwd=str(BUILD_DIR),
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("[ERROR] make falló:\n", r.stderr)
        sys.exit(1)
    print("[OK] Compilación exitosa.")


def ejecutar_simulacion(target, force=False):
    """Ejecuta la simulación para una condición dada."""
    outfile = OUTPUT_DIR / f"data_{target}.txt"
    if outfile.exists() and not force:
        print(f"  [{target}] Datos existentes, cargando...")
        return np.loadtxt(str(outfile))

    exe = BUILD_DIR / "circuitos" / "cgc_ramp"
    print(f"  [{target}] Ejecutando simulación...")
    with open(str(outfile), "w") as f:
        r = subprocess.run([str(exe), target], stdout=f,
                           stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        print(f"  [{target}] ERROR:\n", r.stderr)
        return None
    print(f"  [{target}] OK.")
    return np.loadtxt(str(outfile))


def detect_spikes(time, voltage, threshold=-20.0):
    """Detecta espigas usando find_peaks con umbral de altura."""
    peaks, _ = find_peaks(voltage, height=threshold, distance=10)
    return time[peaks]


def detect_bursts(spike_times, min_spikes=2, max_isi=300.0):
    """
    Se acumulan espigas consecutivas cuyo ISI <= max_isi.
    Cuando el ISI excede max_isi, se cierra la ráfaga actual.
    Solo se conservan ráfagas con >= min_spikes espigas.

    Return: Lista de tuplas (b_start, b_end) para cada burst válido.
    """
    if len(spike_times) < min_spikes:
        return []

    bursts = []
    current_burst_spikes = [spike_times[0]]

    for i in range(1, len(spike_times)):
        isi = spike_times[i] - spike_times[i - 1]
        if isi <= max_isi:
            # Misma ráfaga
            current_burst_spikes.append(spike_times[i])
        else:
            # Fin de ráfaga: evaluar y guardar si cumple min_spikes
            if len(current_burst_spikes) >= min_spikes:
                bursts.append((current_burst_spikes[0], current_burst_spikes[-1]))
            # Iniciar nueva ráfaga
            current_burst_spikes = [spike_times[i]]

    # Última ráfaga pendiente
    if len(current_burst_spikes) >= min_spikes:
        bursts.append((current_burst_spikes[0], current_burst_spikes[-1]))

    return bursts



def compute_intervals(bursts_n1m, bursts_n2v, bursts_n3t):
    """
    Calcula los 15 intervalos de la secuencia CPG + el período ciclo a ciclo.
    
    Definición de un ciclo i (referenciado a N1M):
      - N1_start, N1_end  = bursts_n1m[i]
      - N2_start, N2_end  = bursts_n2v[i]
      - N3_start, N3_end  = bursts_n3t[i]
      - N1_next           = bursts_n1m[i+1][0]  (inicio del siguiente ciclo)
    
    Los intervalos calculados son:
      period:      N1_next - N1_start
      
      Burst Durations (BD):
        N1M_BD:    N1_end - N1_start
        N2v_BD:    N2_end - N2_start
        N3t_BD:    N3_end - N3_start
      
      Intervals:
        N1_N2_int: N2_start - N1_start
        N1_N3_int: N3_start - N1_start
        N2_N3_int: N3_start - N2_start
        N2_N1_int: N1_next  - N2_start
        N3_N1_int: N1_next  - N3_start
        N3_N2_int: N2_start - N3_start
      
      Delays:
        N1_N2_del: N2_start - N1_end
        N1_N3_del: N3_start - N1_end
        N2_N3_del: N3_start - N2_end
        N2_N1_del: N1_next  - N2_end
        N3_N1_del: N1_next  - N3_end
        N3_N2_del: N2_start - N3_end
    
    Retorna un pandas.DataFrame (filas=ciclos, columnas=intervalos),
    o None si no hay suficientes ciclos.
    """
    records = []
    for i in range(len(bursts_n1m) - 1):
        n1_on, n1_off = bursts_n1m[i]
        n1_next = bursts_n1m[i + 1][0]
        period = n1_next - n1_on

        # Encontrar los primeros bursts de N2v y N3t dentro del ciclo actual
        n2_cands = [b for b in bursts_n2v if b[0] >= n1_on and b[0] < n1_next]
        n3_cands = [b for b in bursts_n3t if b[0] >= n1_on and b[0] < n1_next]

        if not n2_cands or not n3_cands:
            continue

        n2_on, n2_off = n2_cands[0]
        n3_on, n3_off = n3_cands[0]

        if period <= 0 or period > 20000:
            continue

        records.append({
            'period':     period,
            # Burst Durations
            'N1M_BD':     n1_off - n1_on,
            'N2v_BD':     n2_off - n2_on,
            'N3t_BD':     n3_off - n3_on,
            # Intervals (onset --> onset)
            'N1_N2_int':  n2_on - n1_on,
            'N1_N3_int':  n3_on - n1_on,
            'N2_N3_int':  n3_on - n2_on,
            'N2_N1_int':  n1_next - n2_on,
            'N3_N1_int':  n1_next - n3_on,
            'N3_N2_int':  n2_on - n3_on,
            # Delays (offset --> onset)
            'N1_N2_del':  n2_on - n1_off,
            'N1_N3_del':  n3_on - n1_off,
            'N2_N3_del':  n3_on - n2_off,
            'N2_N1_del':  n1_next - n2_off,
            'N3_N1_del':  n1_next - n3_off,
            'N3_N2_del':  n2_on - n3_off,
        })

    if len(records) < 3:
        return None

    return pd.DataFrame(records)


def plot_boxplots(df, target):

    df_seconds = df[BOXPLOT_KEYS] / 1000.0
    df_long = df_seconds.melt(var_name='Intervalo', value_name='Duración (s)')
    df_long['Intervalo'] = df_long['Intervalo'].map(LABELS_PRETTY)

    # Orden de aparición
    order = [LABELS_PRETTY[k] for k in BOXPLOT_KEYS]

    # Asignar grupo para colorear por hue
    group_map = {LABELS_PRETTY['period']: 'Period'}
    for k in ['N1M_BD', 'N2v_BD', 'N3t_BD']:
        group_map[LABELS_PRETTY[k]] = 'Burst Duration'
    for k in ['N1_N2_int', 'N1_N3_int', 'N2_N3_int', 'N2_N1_int', 'N3_N1_int', 'N3_N2_int']:
        group_map[LABELS_PRETTY[k]] = 'Interval'
    for k in ['N1_N2_del', 'N1_N3_del', 'N2_N3_del', 'N2_N1_del', 'N3_N1_del', 'N3_N2_del']:
        group_map[LABELS_PRETTY[k]] = 'Delay'
    df_long['Grupo'] = df_long['Intervalo'].map(group_map)

    palette_hue = {
        'Period': 'orange',
        'Burst Duration': 'blue',
        'Interval': 'green',
        'Delay': 'red',
    }

    fig, ax = plt.subplots(figsize=(16, 6))
    sns.boxplot(
        data=df_long,
        x='Intervalo',
        y='Duración (s)',
        hue='Grupo',
        order=order,
        palette=palette_hue,
        showfliers=False,
        linewidth=0.8,
        dodge=False,
        width=0.4,
        legend=True,
        ax=ax,
    )

    ax.set_title(
        f'Variabilidad de intervalos -- Estimulación {target} (CPG+CGC)',
        fontweight='bold', fontsize=11
    )
    ax.set_xlabel('')
    ax.axhline(0, color='gray', lw=0.5, ls='--')
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.tight_layout()

    path = OUTPUT_DIR / f'boxplots_{target}.png'
    fig.savefig(str(path), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  --> {path.name}")


def plot_correlations_pairplot(df, target):
    """
    Grid 5x3 de correlación intervalo vs period con regresión (Fig. 6/8/10).
    Usa sns.regplot en cada celda de un grid manual de 5 filas x 3 columnas.
    Calcula R con scipy.stats.linregress y lo anota en cada subgráfica.
    """
    # 5 filas x 3 columnas: BD(3), Intervals(6), Delays(6) = 15
    rows_keys = [
        ['N1M_BD',     'N2v_BD',     'N3t_BD'],
        ['N1_N2_int',  'N1_N3_int',  'N2_N3_int'],
        ['N2_N1_int',  'N3_N1_int',  'N3_N2_int'],
        ['N1_N2_del',  'N1_N3_del',  'N2_N3_del'],
        ['N2_N1_del',  'N3_N1_del',  'N3_N2_del'],
    ]
    row_titles = ['Burst Duration', 'Interval (onset)', 'Interval (onset)',
                  'Delay (offset)', 'Delay (offset)']

    fig, axes = plt.subplots(5, 3, figsize=(13, 17))
    period = df['period'].values / 1000.0

    for i, (rk, rtitle) in enumerate(zip(rows_keys, row_titles)):
        for j, key in enumerate(rk):
            ax = axes[i, j]
            y = df[key].values / 1000.0

            slope, intercept, r_val, _, _ = linregress(period, y)
            r2 = r_val ** 2

            # Color según R
            if r2 > 0.7:
                color = 'red'
            elif r2 > 0.4:
                color = 'orange'
            else:
                color = 'blue'

            sns.regplot(
                x=period, y=y, ax=ax, ci=None,
                scatter_kws={'s': 12, 'alpha': 0.5, 'color': color, 'edgecolors': 'none'},
                line_kws={'color': 'black', 'linewidth': 1.2},
            )

            lbl = LABELS_PRETTY[key]
            ax.set_title(f'{lbl}  (R²={r2:.3f}, m={slope:.3f})', fontsize=8)
            ax.tick_params(labelsize=7)

            if i == 4:
                ax.set_xlabel('Period (s)', fontsize=8)
            else:
                ax.set_xlabel('')
            if j == 0:
                ax.set_ylabel(f'{rtitle} (s)', fontsize=8)
            else:
                ax.set_ylabel('')

            ax.annotate(
                f'R²={r2:.3f}\nm={slope:.3f}',
                xy=(0.05, 0.95), xycoords='axes fraction',
                fontsize=7, fontweight='bold', va='top', ha='left',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='gray', alpha=0.85),
            )

    fig.suptitle(
        f'Invariantes Dinamicos -- Estimulacion {target} (CPG+CGC)',
        fontsize=12, fontweight='bold', y=0.995
    )
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    path = OUTPUT_DIR / f'correlations_{target}.png'
    fig.savefig(str(path), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  --> {path.name}")


def print_summary(df, target):
    """Imprime tabla resumen de R para una condición."""
    period = df['period'].values
    n_cycles = len(period)
    
    print(f"\n{'=' * 65}")
    print(f"  RESUMEN -- Estimulación {target}  ({n_cycles} ciclos)")
    print()
    print(f"  {'Intervalo':<18} {'Pendiente':>10} {'R²':>8} {'Tipo':>16}")
    print(f"\n{'=' * 65}")

    for key in INTERVAL_KEYS:
        y = df[key].values
        slope, _, r_val, _, _ = linregress(period, y)
        r2 = r_val ** 2
        if r2 > 0.7:
            tipo = " INVARIANTE"
        elif r2 > 0.4:
            tipo = " Moderado"
        else:
            tipo = " Débil"
        print(f"  {LABELS_PRETTY[key]:<18} {slope:>10.4f} {r2:>8.4f} {tipo:>16}")
    print(f"{'=' * 65}")


def plot_comparison(all_results):
    """Figura comparativa: R de cada intervalo para las 4 condiciones."""
    labels = [LABELS_PRETTY[k] for k in INTERVAL_KEYS]

    fig, ax = plt.subplots(figsize=(14, 5))
    x = np.arange(len(INTERVAL_KEYS))
    n_targets = len(all_results)
    width = 0.8 / max(n_targets, 1)

    for idx, (target, df) in enumerate(all_results.items()):
        period = df['period'].values
        r2_vals = []
        for key in INTERVAL_KEYS:
            _, _, r_val, _, _ = linregress(period, df[key].values)
            r2_vals.append(r_val ** 2)
        offset = (idx - n_targets / 2 + 0.5) * width
        ax.bar(
            x + offset, r2_vals, width,
            label=target,
            alpha=0.85, edgecolor='white'
        )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('R²')
    ax.set_title(
        'Comparación de invariantes dinámicos entre condiciones (CPG+CGC)',
        fontweight='bold'
    )
    ax.axhline(0.7, color='red', ls='--', lw=0.8, label='Umbral invariante (R²=0.7)')
    ax.legend(fontsize=8, loc='upper right')
    ax.set_ylim(0, 1.05)
    plt.tight_layout()

    path = OUTPUT_DIR / 'comparacion_R2.png'
    fig.savefig(str(path), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  --> Comparativa: {path.name}")


def analizar_condicion(data, target):
    """Analiza una condición: detecta bursts, calcula intervalos, genera figuras."""
    t = data[:, COL_T]

    mask = t > TRANSIENT_TIME
    ds = data[mask]
    ts = ds[:, COL_T]

    sp_n1m = detect_spikes(ts, ds[:, COL_N1M], SPIKE_THRESHOLDS['N1M'])
    sp_n2v = detect_spikes(ts, ds[:, COL_N2V], SPIKE_THRESHOLDS['N2v'])
    sp_n3t = detect_spikes(ts, ds[:, COL_N3T], SPIKE_THRESHOLDS['N3t'])

    print(f"  Espigas: N1M={len(sp_n1m)}, N2v={len(sp_n2v)}, N3t={len(sp_n3t)}")

    b_n1m = detect_bursts(sp_n1m, MIN_SPIKES['N1M'], MAX_ISI['N1M'])
    b_n2v = detect_bursts(sp_n2v, MIN_SPIKES['N2v'], MAX_ISI['N2v'])
    b_n3t = detect_bursts(sp_n3t, MIN_SPIKES['N3t'], MAX_ISI['N3t'])

    print(f"  Bursts:  N1M={len(b_n1m)}, N2v={len(b_n2v)}, N3t={len(b_n3t)}")

    df = compute_intervals(b_n1m, b_n2v, b_n3t)
    if df is None or len(df) < 5:
        print(f"  [WARN] Insuficientes ciclos para {target}. Saltando.")
        return None

    print(f"  Ciclos válidos: {len(df)}")

    csv_path = OUTPUT_DIR / f'intervalos_{target}.csv'
    df.to_csv(str(csv_path), index=False, float_format='%.4f')
    print(f"  --> {csv_path.name}")

    plot_boxplots(df, target)
    plot_correlations_pairplot(df, target)
    print_summary(df, target)

    return df


def print_ramp_config(target):
    """Imprime tabla de configuración de rampa y corrientes para una condición."""
    cfg = RAMP_CONFIG[target]
    print(f"  Configuración ({cfg['nota']}):")
    print(f"    Corrientes inyectadas (nA):")
    for key in ['I_SO', 'I_N1M', 'I_N2v', 'I_N3t', 'I_CGC']:
        val = cfg[key]
        if val == 'RAMPA':
            print(f"      {key:>6} = RAMPA [{cfg['c_min']} --> {cfg['c_max']}] "
                  f"({cfg['n_steps']} escalones, intervalo={RAMP_INTERVAL} ms)")
        else:
            print(f"      {key:>6} = {val}")
    total_time = 2 * cfg['n_steps'] * RAMP_INTERVAL / 1000.0
    print(f"    Duracion total rampa: {total_time:.1f} s (1 triangulo)")


def main():
    force = '--rerun' in sys.argv

    print(" ESTUDIO DE INVARIANTES DINÁMICOS -- CPG + CGC")
    print(" 4 condiciones: N1M, N3t, SO, CGC")

    compilar()

    all_results = {}
    for target in TARGETS:
        print(f"\n{'-' * 65}")
        print(f" Condición: Rampa en {target}")
        print(f"{'-' * 65}")

        data = ejecutar_simulacion(target, force=force)
        if data is None:
            continue

        print(f"  Datos: {data.shape[0]} muestras, {data.shape[1]} columnas")
        print_ramp_config(target)
        result = analizar_condicion(data, target)
        if result is not None:
            all_results[target] = result

    if len(all_results) > 1:
        print("aaaa")
        plot_comparison(all_results)

    print(f"\n{'=' * 65}")
    print(f" COMPLETADO. Figuras en: {OUTPUT_DIR}")
    print(f"{'=' * 65}")


if __name__ == '__main__':
    main()
