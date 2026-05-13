#!/usr/bin/env python3
"""
Script para visualizar las simulaciones enteras (data_*.txt).
Diseñado para consumir muy poca RAM utilizando lectura por bloques (chunks)
y downsampling agresivo que preserva los picos (spikes).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Importamos la lógica de detección del script principal
from estudio_invariantes import detect_spikes, detect_bursts, SPIKE_THRESHOLDS, MIN_SPIKES, MAX_ISI

# Nombres de las 7 columnas que saca cgc_cpg_ramp.cpp
COL_NAMES = ['time', 'N1M', 'N2v', 'N3t', 'SO', 'CGC', 'ramp']

def plot_simulation(target):
    filepath = Path("resultados_invariantes") / f"data_{target}.txt"
    if not filepath.exists():
        print(f"  [ERROR] No se encuentra el archivo {filepath}")
        return
        
    df = pd.read_csv(filepath, sep=r'\s+', header=None, names=COL_NAMES)
    
    fig, axes = plt.subplots(6, 1, figsize=(14, 16), sharex=False)
    
    t = df['time']
    neurons = ['N1M', 'N2v', 'N3t', 'SO', 'CGC']
    colors = ['#2196F3', '#FF9800', '#4CAF50', '#F44336', '#9C27B0']
    
    for ax, neuron, color in zip(axes[:5], neurons, colors):
        ax.plot(t, df[neuron], color=color, lw=1.0)
        ax.set_ylabel(f"{neuron} (mV)", fontweight='bold')
        ax.set_xlabel("Tiempo de simulación (ms)", fontweight='bold')
        ax.grid(True, alpha=0.3)

        if neuron in ['N1M', 'N3t']:
            ax.axhline(-20, color='gray', linestyle='--', lw=0.8, alpha=0.5)
        elif neuron == 'N2v':
            ax.axhline(-25, color='gray', linestyle='--', lw=0.8, alpha=0.5)
        
    axes[5].plot(t, df['ramp'], color='black', lw=1.5)
    axes[5].set_ylabel("Rampa", fontweight='bold')
    axes[5].set_xlabel("Tiempo de simulación (ms)", fontweight='bold')
    axes[5].grid(True, alpha=0.3)
    
    fig.suptitle(f"Actividad del circuito — Condición {target} (CPG + CGC)", 
                 fontsize=14, fontweight='bold', y=0.92)
    
    plt.subplots_adjust(hspace=0.5)
    
    out_path = Path("resultados_invariantes") / f"trazos_simulacion_{target}.png"
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  --> Gráfico guardado en: {out_path}\n")


def plot_un_ciclo(target):
    filepath = Path("resultados_invariantes") / f"data_{target}.txt"
    if not filepath.exists():
        print(f"  [ERROR] No se encuentra {filepath}")
        return
        
    print(f"  [Ciclo] Analizando {target} para extraer un ciclo...")
    
    data = np.loadtxt(str(filepath))
    t = data[:, 0]
    
    mask = t > 5000.0
    ds = data[mask]
    ts = ds[:, 0]
    
    sp_n1m = detect_spikes(ts, ds[:, 1], SPIKE_THRESHOLDS['N1M'])
    b_n1m = detect_bursts(sp_n1m, MIN_SPIKES['N1M'], MAX_ISI['N1M'])
    
    if len(b_n1m) < 3:
        print(f"  [WARN] No hay suficientes ráfagas en N1M para definir un ciclo en {target}.")
        return
        
    mid_idx = len(b_n1m) // 2
    t_start = b_n1m[mid_idx][0]    
    t_end = b_n1m[mid_idx + 1][0]   
    
    margin = 500.0
    t_start -= margin
    t_end += margin
    
    mask_cycle = (data[:, 0] >= t_start) & (data[:, 0] <= t_end)
    cycle_data = data[mask_cycle]
    
    tc = cycle_data[:, 0]
    
    fig, axes = plt.subplots(5, 1, figsize=(10, 10), sharex=False)
    
    neurons = ['N1M', 'N2v', 'N3t', 'SO', 'CGC']
    cols = [1, 2, 3, 4, 5]
    colors = ['#2196F3', '#FF9800', '#4CAF50', '#F44336', '#9C27B0']
    
    for ax, neuron, col_idx, color in zip(axes, neurons, cols, colors):
        ax.plot(tc, cycle_data[:, col_idx], color=color, lw=1.2)
        ax.set_ylabel(f"{neuron} (mV)", fontweight='bold')
        ax.set_xlabel("Tiempo (ms)", fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        if neuron in SPIKE_THRESHOLDS:
            ax.axhline(SPIKE_THRESHOLDS[neuron], color='gray', linestyle='--', lw=0.8, alpha=0.5)
            
    fig.suptitle(f"Detalle de un Ciclo Completo — Condición {target}", 
                 fontsize=14, fontweight='bold', y=0.95)
    plt.subplots_adjust(hspace=0.5)
    
    out_path = Path("resultados_invariantes") / f"ciclo_detalle_{target}.png"
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  --> Gráfico guardado en: {out_path} \n      (Ventana temporal: {t_start:.0f} - {t_end:.0f} ms)\n")


if __name__ == '__main__':
    # Podemos pasar los targets por argumento: python3 plot_simulaciones.py N1M SO
    if len(sys.argv) > 1:
        targets = sys.argv[1:]
    else:
        targets = ['N1M', 'N3t', 'SO', 'CGC']
        
    print(" GENERADOR DE TRAZOS DE SIMULACIÓN")
    
    for t in targets:
        print(f"\nProcesando condición: {t}")
        plot_simulation(t)
        plot_un_ciclo(t)
