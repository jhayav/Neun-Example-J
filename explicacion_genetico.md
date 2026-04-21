# Explicación Detallada: Algoritmo Genético para la Parametrización de la CGC

Este documento detalla el diseño, la justificación biológica y la implementación del algoritmo genético (GA) utilizado para encontrar los parámetros sinápticos óptimos de la neurona CGC (Cerebral Giant Cell) sobre el circuito CPG (Central Pattern Generator) del motor trifásico de *Lymnaea stagnalis*.

---

## 1. Motivación del Enfoque Híbrido

Las conexiones sinápticas de las neuronas del CPG (N1M, N2v, N3t, SO) han sido profundamente estudiadas y parametrizadas en estudios previos (Vavoulis et al., 2007). Sin embargo, los pesos sinápticos exactos (conductancias `Gsyn`) y los estímulos de la CGC para lograr un comportamiento fidedigno no están modelados matemáticamente con el mismo rigor; sólo disponemos de sus **efectos fenomenológicos** observados en laboratorio.

Para cerrar esta brecha, implementamos un enfoque de **Ajuste de Parámetros Guiado por Datos (Parameter Fitting)** usando un Algoritmo Genético (AG). El simulador, sumamente rápido pero rígido, está escrito en **C++** (`circuito_gen.cpp`). El proceso de optimización, que requiere un análisis matemático avanzado de picos (spikes), métricas de regularidad y lógica de selección biológica, está orquestado en **Python** (`genetic_search_bio.py`).

---

## 2. Arquitectura del Sistema Evolutivo

### 2.1. El Genotipo (El Cromosoma)
Cada individuo en nuestra población representa un conjunto posible de configuraciones de la red. Nuestro "cromosoma" consta de 5 "genes" (valores flotantes), confinados en límites con sentido físico:

1. **`Gsyn_CGC_N1M`** [0.0 - 50.0 mS/cm²]: Fuerza de excitación directa sobre N1M.
2. **`Gsyn_CGC_N2v`** [0.0 - 50.0 mS/cm²]: Fuerza de excitación sobre N2v.
3. **`Gsyn_CGC_SO`** [0.0 - 50.0 mS/cm²]: Fuerza de excitación sobre SO.
4. **`Gsyn_CGC_N3t`** [0.0 - 50.0 mS/cm²]: Fuerza de excitación sobre N3t.
5. **`I_ext_CGC`** [-4.0 - 4.0 nA]: Corriente basal inyectada en la CGC para emular su drive tónico.

### 2.2. Parámetros del Algoritmo
- **Población (N=30):** Se evalúan 30 configuraciones de red distintas en cada generación.
- **Generaciones (Gen=100):** Tiempo suficiente para asegurar convergencia profunda sin eternizar la ejecución.
- **Elitismo (Top 5%):** El mejor individuo de cada generación sobrevive intacto a la siguiente, garantizando que el AG nunca "olvide" una solución excelente encontrada previamente. En una población de 30, esto protege al *Top 1*.
- **Selección (Torneo Binario):** Se eligen dos individuos al azar y el que tiene mejor *fitness* se convierte en "padre". Esto mantiene la presión evolutiva mientras preserva la diversidad.
- **Cruce (Crossover simple 70%):** Con un 70% de probabilidad, dos padres mezclan sus genes en un punto aleatorio, produciendo dos configuraciones híbridas.
- **Mutación (Mutación Gaussiana 20%):** Con un 20% de probabilidad, un gen sufre una perturbación basada en una distribución normal (desviación fina de 0.15). Esto es fundamental para hacer ajustes precisos y evitar mínimos locales.

---

## 3. El Corazón del Sistema: La Función de Fitness

La función de evaluación es el núcleo donde las **restricciones neurobiológicas reales** descritas en la literatura se traducen en ecuaciones de recompensa (bonus) y castigo (penalización). 

La función toma como *"input"* las trazas de voltaje devueltas por C++ a lo largo de 10 segundos (10,000 ms), y genera un *Score* numérico global.

### Regla Biológica 1: Firing Rate Fisiológico de la CGC
**Biología:** *"During feeding the CGCs fire maximally in the 7 to 20 spike/minute range... above the threshold level of firing, the CGCs also influence the frequency."*
* **Implementación:**
  Si medimos en spikes por minuto (SPM), en 10,000 ms (10 segundos, que es un sexto de minuto) deberíamos observar matemáticamente entre $1.16$ y $3.33$ disparos.
  * **Castigo extremo (-2000):** Si hay $0$ picos, la CGC está apagada previendo la ingesta (quiescence / locomotion). Si hay $>5$, su frecuencia es antifisiológica.
  * **Recompensa (+500):** Si está perfectamente en el rango de 1 a 4 picos de voltaje por encima de -20mV en la ventana de 10 segundos.

### Regla Biológica 2: Disparo Tónico Regular
**Biología:** *"Continuous or tonic spiking activity in the CGCs provides a background of excitatory modulation."*
* **Implementación:**
  No buscamos ráfagas (*bursting*), sino un disparo regular de fondo.
  Para ello medimos el **ISI** (Inter-Spike Interval). Se calcula la desviación estándar del ISI dividida por su media (Coeficiente de Variación, CV). Cuanto inferior sea la varianza, más regular/constante es el ritmo de la CGC.
  * **Bono:** `score += 200 * (1.0 - CV)`. Penaliza cronometrías caóticas y premia cadencias exactas.

### Regla Biológica 3: El Papel Modulador ("Gating") en N1M y N2v
**Biología:** *"The CGCs reduce the threshold for plateauing in both of these neuron types... making them more likely to respond... This is usually too weak to initiate plateaus [on its own in N2v] but if... depolarized... plateaus are initiated as is the case with N1Ms."*
* **Implementación:**
  La CGC no *crea* el patrón por sí misma, sino que empuja a `N1M` y `N2v` lo suficientemente alto en voltaje para que sus dinámicas intrínsecas hagan "plateau".
  * **N1M:** Requerimos que rompa el umbral de $-40 mV$. Si el voltaje máximo se queda corto se penaliza brutalmente (`-1000`). Si lo supera, se suma un bono escalado `(Vmax - (-40.0)) * 10`.
  * **N2v:** Igual, pero con un umbral ligeramente más fácil ($-45 mV$), penalizando con `-500` si ni siquiera se despolariza durante la simulación. Recompensa lineal por llegar a valores altos de plateau.

### Regla Biológica 4: Comprobación de Integridad del Circuito
**Biología:** Toda intervención no debe suprimir ni corromper las propiedades oscilatorias del CPG trifásico.
* **Implementación:** Contamos de forma agresiva mediante la función `scipy.signal.find_peaks` que `N1M` y `N3t` presenten al menos dos ráfagas completas en sus ritmos de oscilación lenta.
  * **Castigo por rotura (-1500):** Si las conductancias introducidas por CGC son tan erráticas o fuertemente inhibitorias (por dinámicas indirectas de red) que colapsan el ritmo y reducen todo a estasis o ruido, el ADN se rechaza.

---

## 4. Resumen

El Algoritmo Genético opera como si se estuvieran haciendo miles de parches electrofisiológicos virtuales. 

El modelo desecha instintivamente los pesos siápticos que inhiben el motor básico o que empujan a la neurona CGC a rangos antifisiológicos (como actuar como oscilador en lugar de neuromodulador). Generación tras generación, los genes irán convergiendo hasta encontrar exactamente **ese punto dulce (sweet spot)** en el cual una inyección corriente muy tenue sobre la CGC propicia un tren constante de picos, el cual empuja delicadamente los voltajes basales de las CPGs permitiéndoles hacer bursting.
