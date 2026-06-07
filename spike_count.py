import numpy as np
import sys

FILE        = 'cgc_1.txt'   # archivo de datos
V_COL       = 9             # columna del voltaje CGC (0-indexed en el archivo raw)
THRESHOLD   = 0.0           # mV — umbral de detección de spike
REFRACTORY  = 5.0           # ms — periodo refractario mínimo entre spikes

if len(sys.argv) > 1:
    FILE = sys.argv[1]

print(f"Leyendo '{FILE}', columna V_CGC = {V_COL}...")

data = np.loadtxt(FILE, usecols=[0, V_COL])
t_ms = data[:, 0]   # tiempo en ms
v    = data[:, 1]   # voltaje CGC en mV

# Detección de spikes: cruce ascendente del umbral con periodo refractario
above      = v >= THRESHOLD
crossings  = np.where(~above[:-1] & above[1:])[0] + 1

spike_times = []
last_spike  = -np.inf
for idx in crossings:
    t = t_ms[idx]
    if t - last_spike >= REFRACTORY:
        spike_times.append(t)
        last_spike = t

spike_times = np.array(spike_times)
n_spikes    = len(spike_times)

t_total_min = (t_ms[-1] - t_ms[0]) / 60000.0   # ms -> min

print(f"\nDuración total   : {t_ms[-1]/1000:.1f} s  ({t_total_min:.3f} min)")
print(f"Spikes detectados: {n_spikes}")
print(f"Tasa global      : {n_spikes / t_total_min:.1f} spikes/min")

# Tasa por ventanas de 1 minuto (si la simulación es suficientemente larga)
window_ms = 60_000.0
t_start   = t_ms[0]
t_end     = t_ms[-1]
n_windows = int((t_end - t_start) / window_ms)

if n_windows >= 2:
    print(f"\nTasa por ventanas de 1 min:")
    for i in range(n_windows):
        w0 = t_start + i * window_ms
        w1 = w0 + window_ms
        n  = np.sum((spike_times >= w0) & (spike_times < w1))
        print(f"  [{w0/1000:.0f}s – {w1/1000:.0f}s]  →  {n} spikes/min")
elif n_spikes > 0:
    # Ventanas de 10 s para simulaciones cortas
    window_ms = 10_000.0
    n_windows = int((t_end - t_start) / window_ms)
    if n_windows >= 2:
        print(f"\nTasa por ventanas de 10 s:")
        for i in range(n_windows):
            w0  = t_start + i * window_ms
            w1  = w0 + window_ms
            n   = np.sum((spike_times >= w0) & (spike_times < w1))
            rate = n / (window_ms / 60000.0)
            print(f"  [{w0/1000:.0f}s – {w1/1000:.0f}s]  →  {rate:.1f} spikes/min")

if n_spikes > 1:
    isis = np.diff(spike_times)          # inter-spike intervals (ms)
    print(f"\nISI medio : {isis.mean():.1f} ms  →  {60000/isis.mean():.1f} spikes/min instantáneo")
    print(f"ISI mín   : {isis.min():.1f} ms")
    print(f"ISI máx   : {isis.max():.1f} ms")
