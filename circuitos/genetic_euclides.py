"""
genetic_euclides.py - Algoritmo Genético con DEAP

Optimiza los 12 parámetros de las sinapsis CGC->CPG (gsyn, esyn, tau_syn × 4)
minimizando la distancia euclidiana entre la salida de circuito_gen y la
referencia de cpg_completo (CPG aislado sin CGC).

Neuronas comparadas: N1M, N2v, N3t, SO (voltaje somático)

Librería: DEAP (Distributed Evolutionary Algorithms in Python)
https://deap.readthedocs.io/
"""

import subprocess
import numpy as np
import os
import random
import sys

from deap import base, creator, tools, algorithms

# CONFIGURACIÓN
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GEN_EXECUTABLE  = os.path.join(SCRIPT_DIR, "../build/circuitos/circuito_gen")
REF_EXECUTABLE  = os.path.join(SCRIPT_DIR, "../build/circuitos/CPG")

POPULATION_SIZE = 30
GENERATIONS = 100
MUTATION_RATE = 0.2    # Probabilidad de mutar cada gen (indpb en DEAP)
CROSSOVER_RATE = 0.7   # Probabilidad de cruce (cxpb en DEAP)
MUTATION_SIGMA_FRAC = 0.05  # Sigma de mutación como fracción del rango

# Subsampling: cada SUBSAMPLE_STEP pasos (paso de simulación = 0.01ms)
# 100 -> paso efectivo de 1ms, ~10000 puntos de comparación
SUBSAMPLE_STEP = 100

# 12 genes: 3 por cada sinapsis CGC->neurona (gsyn, esyn, tau_syn)
# Todas las sinapsis son EXCITATORIAS -> esyn alrededor de 0 mV
GENE_NAMES = [
    "gsyn_CGC_N1M", "esyn_CGC_N1M", "tau_CGC_N1M",
    "gsyn_CGC_N2v", "esyn_CGC_N2v", "tau_CGC_N2v",
    "gsyn_CGC_SO",  "esyn_CGC_SO",  "tau_CGC_SO",
    "gsyn_CGC_N3t", "esyn_CGC_N3t", "tau_CGC_N3t",
]

GENE_BOUNDS = [
    # CGC -> N1M
    (0.0,  10.0),   # gsyn
    (-10.0, 10.0),  # esyn
    (10.0, 500.0),  # tau_syn (ms)
    # CGC -> N2v
    (0.0,  10.0),
    (-10.0, 10.0),
    (10.0, 500.0),
    # CGC -> SO
    (0.0,  10.0),
    (-10.0, 10.0),
    (10.0, 500.0),
    # CGC -> N3t
    (0.0,  10.0),
    (-10.0, 10.0),
    (10.0, 500.0),
]

LOW_BOUNDS  = [b[0] for b in GENE_BOUNDS]
HIGH_BOUNDS = [b[1] for b in GENE_BOUNDS]
N_GENES = len(GENE_BOUNDS)

# Sigmas de mutación gaussiana (proporcional al rango de cada gen)
MUTATION_SIGMAS = [MUTATION_SIGMA_FRAC * (hi - lo) for lo, hi in GENE_BOUNDS]


# SEÑALES DE REFERENCIA (variable global, se carga una sola vez)
REF_SIGNALS = None  # Se inicializa en main()


def run_reference():
    """Ejecuta cpg_completo (CPG) una sola vez y devuelve las señales de referencia.

    Salida de CPG: 23 columnas
    Col 0: tiempo
    Col 1: V_N1M_soma, Col 2: V_N1M_axon
    Col 3: V_N2v_soma, Col 4: V_N2v_axon
    Col 5: V_N3t_soma, Col 6: V_N3t_axon
    Col 7: V_SO_soma,  Col 8: V_SO_axon
    ... (corrientes y gating)

    devuelve: (t, v_n1m, v_n2v, v_n3t, v_so) submuestreados
    """
    print("Ejecutando simulación de referencia (cpg_completo)...")
    try:
        res = subprocess.run([REF_EXECUTABLE], capture_output=True, text=True,
                             check=True, timeout=120)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        print(f"ERROR ejecutando referencia: {e}")
        sys.exit(1)

    lines = res.stdout.strip().split('\n')
    data = []
    for line in lines:
        if line.strip():
            vals = [float(x) for x in line.split()]
            data.append(vals)

    data = np.array(data)

    # Subsampling
    data = data[::SUBSAMPLE_STEP]

    t      = data[:, 0]
    v_n1m  = data[:, 1]   # V_N1M_soma
    v_n2v  = data[:, 3]   # V_N2v_soma
    v_n3t  = data[:, 5]   # V_N3t_soma
    v_so   = data[:, 7]   # V_SO_soma

    print(f"  Referencia cargada: {len(t)} puntos (subsampled)")
    return t, v_n1m, v_n2v, v_n3t, v_so


def run_gen_simulation(genes):
    """Ejecuta circuito_gen con los 12 parámetros.

    Salida de circuito_gen: 6 columnas
    Col 0: tiempo, Col 1: V_CGC, Col 2: V_N1M, Col 3: V_N2v, Col 4: V_N3t, Col 5: V_SO

    devuelve: (v_n1m, v_n2v, v_n3t, v_so) submuestreados, o None si falla
    """
    cmd = [GEN_EXECUTABLE] + [str(g) for g in genes]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True,
                             check=True, timeout=120)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        return None

    lines = res.stdout.strip().split('\n')
    data = []
    for line in lines:
        if line.strip():
            vals = [float(x) for x in line.split()]
            data.append(vals)

    if len(data) == 0:
        return None

    data = np.array(data)

    if np.isnan(data).any() or np.isinf(data).any():
        return None

    # Subsampling
    data = data[::SUBSAMPLE_STEP]

    v_n1m = data[:, 2]
    v_n2v = data[:, 3]
    v_n3t = data[:, 4]
    v_so  = data[:, 5]

    return v_n1m, v_n2v, v_n3t, v_so



def evaluate(individual):
    """Calcula fitness = -distancia euclidiana total entre simulación y referencia.

    DEAP requiere que devuelva una tupla.
    """
    global REF_SIGNALS
    ref_n1m, ref_n2v, ref_n3t, ref_so = REF_SIGNALS

    result = run_gen_simulation(individual)
    if result is None:
        return (-1e9,)  # Penalización máxima

    v_n1m, v_n2v, v_n3t, v_so = result

    # Asegurar misma longitud
    n = min(len(v_n1m), len(ref_n1m))

    # Distancia euclidiana por neurona
    dist_n1m = np.sqrt(np.sum((ref_n1m[:n] - v_n1m[:n])**2))
    dist_n2v = np.sqrt(np.sum((ref_n2v[:n] - v_n2v[:n])**2))
    dist_n3t = np.sqrt(np.sum((ref_n3t[:n] - v_n3t[:n])**2))
    dist_so  = np.sqrt(np.sum((ref_so[:n]  - v_so[:n])**2))

    dist_total = dist_n1m + dist_n2v + dist_n3t + dist_so

    return (-dist_total,)


#  OPERADOR DE MUTACIÓN PERSONALIZADO

def mutate_bounded(individual, mu, sigma, indpb):
    """Mutación gaussiana con clipping a los límites de cada gen."""
    for i in range(len(individual)):
        if random.random() < indpb:
            individual[i] += random.gauss(mu, sigma[i])
            individual[i] = max(LOW_BOUNDS[i], min(individual[i], HIGH_BOUNDS[i]))
    return (individual,)



def setup_deap():
    """Configura el framework DEAP: tipos, operadores, estadísticas."""

    # FitnessMax porque fitness = -distancia (maximizar = minimizar distancia)
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()

    # Generador de genes: un float uniforme por cada gen dentro de sus límites
    def init_gene(i):
        return random.uniform(LOW_BOUNDS[i], HIGH_BOUNDS[i])

    def init_individual():
        return creator.Individual([init_gene(i) for i in range(N_GENES)])

    toolbox.register("individual", init_individual)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # Operadores genéticos
    toolbox.register("evaluate", evaluate)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", mutate_bounded, mu=0, sigma=MUTATION_SIGMAS,
                     indpb=MUTATION_RATE)
    toolbox.register("select", tools.selTournament, tournsize=3)

    # Estadísticas
    stats = tools.Statistics(lambda ind: ind.fitness.values[0])
    stats.register("min_dist", lambda x: -max(x))       # Distancia mínima
    stats.register("avg_dist", lambda x: -np.mean(x))   # Distancia media
    stats.register("max_dist", lambda x: -min(x))       # Distancia máxima

    # Hall of Fame (guarda los N mejores individuos de toda la evolución)
    hof = tools.HallOfFame(5)

    return toolbox, stats, hof



def main():
    global REF_SIGNALS

    print("AG con DEAP - Distancia Euclidiana CPG-CGC")

    # 1. Obtener señales de referencia (cpg_completo)
    t_ref, ref_n1m, ref_n2v, ref_n3t, ref_so = run_reference()
    REF_SIGNALS = (ref_n1m, ref_n2v, ref_n3t, ref_so)

    # 2. Configurar DEAP
    toolbox, stats, hof = setup_deap()

    print(f"\nPoblación: {POPULATION_SIZE} | Generaciones: {GENERATIONS}")
    print(f"Genes: {N_GENES} | Cruce: {CROSSOVER_RATE} | Mutación indpb: {MUTATION_RATE}")
    print()

    # 3. Crear población inicial
    pop = toolbox.population(n=POPULATION_SIZE)

    # 4. Ejecutar algoritmo genético con eaSimple
    #    eaSimple = selección + cruce + mutación + elitismo (vía HallOfFame)
    pop, logbook = algorithms.eaSimple(
        pop, toolbox,
        cxpb=CROSSOVER_RATE,
        mutpb=1.0,           # Probabilidad de que un individuo sea mutado
                             # La prob. por gen la controla indpb en mutate_bounded
        ngen=GENERATIONS,
        stats=stats,
        halloffame=hof,
        verbose=True
    )

    # 5. Resultados finales
    print("RESULTADOS FINALES")

    best = hof[0]
    best_fit = best.fitness.values[0]

    print(f"\nDistancia euclidiana mínima: {-best_fit:.4f}")
    print(f"\nMejores parámetros:")

    output_path = os.path.join(SCRIPT_DIR, "mejores_parametros_euclides.txt")
    with open(output_path, "w") as f:
        f.write(f"# Mejores parámetros - AG DEAP - Distancia euclidiana\n")
        f.write(f"# Fitness (neg dist): {best_fit:.6f}\n")
        f.write(f"# Distancia total: {-best_fit:.6f}\n\n")
        for name, val in zip(GENE_NAMES, best):
            print(f"  {name:20s} = {val:.6f}")
            f.write(f"{name}={val}\n")

    print(f"\nGuardado en: {output_path}")

    # Comando para ejecutar con los mejores parámetros
    cmd = f"{GEN_EXECUTABLE} " + " ".join([f"{x:.6f}" for x in best])
    print(f"\nComando de ejecución:\n  {cmd}")

    # 6. Mostrar top 5 del Hall of Fame
    print(f"\n--- Hall of Fame (top {len(hof)}) ---")
    for rank, ind in enumerate(hof):
        dist = -ind.fitness.values[0]
        print(f"  #{rank+1}: dist={dist:.2f}  genes={[round(g, 4) for g in ind]}")


if __name__ == "__main__":
    main()
