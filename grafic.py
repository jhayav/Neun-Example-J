import numpy as np
import matplotlib.pyplot as plt
import time

print("Cargando 'sim.txt' de forma optimizada para archivos grandes...")
start_time = time.time()

# 0: tiempo
# 1, 3, 5, 7, 9: V_N1M_s, V_N2v_s, V_N3t_s, V_SO_s, V_CGC
# 10, 11, 12, 13, 14: Corrientes N1M-N2v, N2v-N1M, N1M-N3t, N3t-N1M, N2v-N3t
# 15, 16, 17: Corrientes N2v-SO, SO-N1M, SO-N2v
# 18, 19, 20, 21: Corrientes CGC-N1M, CGC-N2v, CGC-SO, CGC-N3t
# 22, 23, 24, 25, 26, 28: p_N1M, p_N2v, q_N2v, p_N3t, q_N3t, h_CGC
usecols = [0, 1, 3, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 28]

stride = 100

try:
    import pandas as pd
    df = pd.read_csv('sim.txt', sep=r'\s+', header=None, usecols=usecols, 
                     skiprows=lambda x: x % stride != 0)
    data = df.values
except ImportError:
    print("Pandas no encontrado, usando lectura línea a línea (puede tomar unos minutos)...")
    parsed_data = []
    with open('sim.txt', 'r') as f:
        for i, line in enumerate(f):
            if i % stride == 0:
                parts = line.split()
                if len(parts) >= 29:
                    row = [float(parts[j]) for j in usecols]
                    parsed_data.append(row)
    data = np.array(parsed_data)

print(f"Carga completada en {time.time() - start_time:.2f} segundos. Puntos a graficar: {data.shape[0]}")

t_cpg = data[:, 0]

# Voltajes somáticos
v_n1m_soma_cpg = data[:, 1]
v_n2v_soma_cpg = data[:, 2]
v_n3t_soma_cpg = data[:, 3]
v_so_soma_cpg  = data[:, 4]
v_cgc_soma_cpg = data[:, 5]


t_max_plot_cpg = t_cpg[-1]
mask_cpg = (t_cpg >= 0) & (t_cpg <= t_max_plot_cpg)

fig, axes = plt.subplots(5, 1, figsize=(15, 13), sharex=True)

# 1. CGC (Cerebral Giant Cell - Neurona Modulatoria)
ax0 = axes[0]
ax0.plot(t_cpg[mask_cpg]/1000, v_cgc_soma_cpg[mask_cpg], 'C1', linewidth=0.8, label=r'$V_S$ CGC')
ax0.axhline(y=-60, color='gray', linestyle='--', alpha=0.5, label=r'$E_{leak}$')
ax0.set_ylabel(r'$V_{CGC}$ (mV)')
ax0.set_title('CGC: Cerebral Giant Cell (Neurona Modulatoria)', fontsize=12, fontweight='bold')
ax0.legend(loc='upper right', fontsize=8)
ax0.set_ylim([-90, 60])

# 2. SO (Modulador de Frecuencia)
ax1 = axes[1]
ax1.plot(t_cpg[mask_cpg]/1000, v_so_soma_cpg[mask_cpg], 'purple', linewidth=0.8, label=r'$V_S$ SO')
ax1.axhline(y=-67, color='gray', linestyle='--', alpha=0.5, label=r'$E_{leak}$')
ax1.set_ylabel(r'$V_{SO}$ (mV)')
ax1.set_title('SO: Slow Oscillator (Modulador de Frecuencia)', fontsize=12, fontweight='bold')
ax1.legend(loc='upper right', fontsize=8)
ax1.set_ylim([-90, 20])

# 3. N1M (Fase 1 - Protracción)
ax2 = axes[2]
ax2.plot(t_cpg[mask_cpg]/1000, v_n1m_soma_cpg[mask_cpg], 'b', linewidth=0.8, label=r'$V_S$ N1M')
ax2.axhline(y=-67, color='gray', linestyle='--', alpha=0.5, label=r'$E_{leak}$')
ax2.axhline(y=-30, color='salmon', linestyle=':', alpha=0.7, label=r'$E_{ACh}$')
ax2.set_ylabel(r'$V_{N1M}$ (mV)')
ax2.set_title('N1M: Interneurona de Protracción (Fase P)', fontsize=12, fontweight='bold')
ax2.legend(loc='upper right', fontsize=8)
ax2.set_ylim([-90, 20])

# 4. N2v (Fase 2 - Rasp)
ax3 = axes[3]
ax3.plot(t_cpg[mask_cpg]/1000, v_n2v_soma_cpg[mask_cpg], 'r', linewidth=0.8, label=r'$V_S$ N2v')
ax3.axhline(y=-67, color='gray', linestyle='--', alpha=0.5, label=r'$E_{leak}$')
ax3.axhline(y=55, color='orange', linestyle=':', alpha=0.7, label=r'$E_{NaL}$')
ax3.set_ylabel(r'$V_{N2v}$ (mV)')
ax3.set_title('N2v: Interneurona de Rasp (Fase R)', fontsize=12, fontweight='bold')
ax3.legend(loc='upper right', fontsize=8)
ax3.set_ylim([-80, 60])

# 5. N3t (Fase 3 - Swallow)
ax4 = axes[4]
ax4.plot(t_cpg[mask_cpg]/1000, v_n3t_soma_cpg[mask_cpg], 'g', linewidth=0.8, label=r'$V_S$ N3t')
ax4.axhline(y=-67, color='gray', linestyle='--', alpha=0.5, label=r'$E_{leak}$')
ax4.axhline(y=80, color='brown', linestyle=':', alpha=0.7, label=r'$E_{Ca}$')
ax4.set_ylabel(r'$V_{N3t}$ (mV)')
ax4.set_title('N3t: Interneurona de Swallow (Fase S)', fontsize=12, fontweight='bold')
ax4.legend(loc='upper right', fontsize=8)
ax4.set_ylim([-90, 100])
ax4.set_xlim([-2, t_max_plot_cpg/1000])

# Marcas de tiempo de la simulación (según cpg_cgc.cpp)
t_stim_start = 100 / 1000   # 0.1 s
t_stim_end = 45000 / 1000   # 45.0 s

for ax in axes:
    ax.axvline(x=t_stim_start, color='k', linestyle='--', alpha=0.4)
    ax.axvline(x=t_stim_end, color='k', linestyle='--', alpha=0.4)

# Etiqueta para la leyenda en el primer subplot
axes[0].axvline(x=t_stim_start, color='k', linestyle='--', alpha=0.4, label='Estimulación')
axes[0].legend(loc='upper right', fontsize=8)

# Configurar marcas explícitas en el eje X cada 5 segundos
ticks = np.arange(0, (t_max_plot_cpg/1000) + 1, 5)

for ax in axes:
    ax.set_xticks(ticks)
    ax.tick_params(labelbottom=True)

ax4.set_xlabel('Tiempo (s)', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig('cpg_cgc_completo.png', dpi=150, bbox_inches='tight')
print("\nGráfica completada y guardada como 'cpg_cgc_completo.png'.")

print(f"\nRangos de voltaje:")
print(f"  CGC soma: {v_cgc_soma_cpg.min():.2f} a {v_cgc_soma_cpg.max():.2f} mV")
print(f"  SO soma: {v_so_soma_cpg.min():.2f} a {v_so_soma_cpg.max():.2f} mV")
print(f"  N1M soma: {v_n1m_soma_cpg.min():.2f} a {v_n1m_soma_cpg.max():.2f} mV")
print(f"  N2v soma: {v_n2v_soma_cpg.min():.2f} a {v_n2v_soma_cpg.max():.2f} mV")
print(f"  N3t soma: {v_n3t_soma_cpg.min():.2f} a {v_n3t_soma_cpg.max():.2f} mV")
