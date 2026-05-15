/*************************************************************
 * cgc_freq_search.cpp - Búsqueda de corriente para frecuencia
 *                       de disparo deseada en la CGC aislada
 *
 * La CGC se dispara en el rango fisiológico de 7–20 veces/min.
 * Este programa busca, mediante búsqueda binaria, qué corriente
 * inyectada produce exactamente 7, 20, 42 y 60 spikes/min.
 *
 * Estrategia:
 *   1. Simular la CGC aislada durante SIM_TIME ms con corriente
 *      constante I_inj. Los primeros TRANSIENT ms se descartan.
 *   2. Contar cruces ascendentes del umbral de spike (V_thr).
 *   3. Extrapolar la tasa a spikes/min.
 *   4. Búsqueda binaria: ajustar I hasta alcanzar la tasa objetivo.
 *
 * Compilación: añadido a CMakeLists.txt (target cgc_freq_search)
 * Uso:          ./cgc_freq_search
 *
 *************************************************************/

#include <DifferentialNeuronWrapper.h>
#include <VavoulisCGCModel.h>
#include <RungeKutta4.h>
#include <SystemWrapper.h>

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

// ─── Tipos ───────────────────────────────────────────────────────────────────
typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<VavoulisCGCModel<double>>, Integrator> CGCNeuron;

// ─── Parámetros de simulación ─────────────────────────────────────────────────
static const double STEP       = 0.01;    // Paso RK4 (ms)
static const double TRANSIENT  = 10000.0; // Transitorio a descartar (ms) = 10 s
static const double SIM_WINDOW = 60000.0; // Ventana de conteo (ms) = 60 s
static const double SIM_TIME   = TRANSIENT + SIM_WINDOW; // Duración total (ms)
static const double V_THR      = -20.0;  // Umbral de detección de spike (mV)

// Número de iteraciones de búsqueda binaria
static const int MAX_ITER = 50;
// Tolerancia: si |tasa_medida - tasa_objetivo| < TOL se considera convergido
static const double TOL_RATE = 0.05;  // spikes/min
// Tolerancia en corriente (para terminar si el intervalo es muy pequeño)
static const double TOL_I = 1e-8;

// ─── Función auxiliar: construir y resetear la neurona CGC ────────────────────
static void init_cgc(CGCNeuron::ConstructorArgs& args) {
    args.params[CGCNeuron::cm]  = 1.0;
    args.params[CGCNeuron::vna] = 55.0;
    args.params[CGCNeuron::vk]  = -90.0;
    args.params[CGCNeuron::vca] = 80.0;

    args.params[CGCNeuron::Gnat] = 1.68;
    args.params[CGCNeuron::Gnap] = 0.44;
    args.params[CGCNeuron::Ga]   = 18.82;
    args.params[CGCNeuron::Gd]   = 1.20;
    args.params[CGCNeuron::Glva] = 0.01;
    args.params[CGCNeuron::Ghva] = 1.03;

    args.params[CGCNeuron::vh_h]    = -56.43;
    args.params[CGCNeuron::vs_h]    = -8.41;
    args.params[CGCNeuron::tau0_h]  = 778.82;
    args.params[CGCNeuron::delta_h] = 0.03;

    args.params[CGCNeuron::vh_r]    = -47.03;
    args.params[CGCNeuron::vs_r]    = 20.55;
    args.params[CGCNeuron::tau0_r]  = 4.01;
    args.params[CGCNeuron::delta_r] = 1.00;

    args.params[CGCNeuron::vh_a]    = -36.37;
    args.params[CGCNeuron::vs_a]    = 8.72;
    args.params[CGCNeuron::tau0_a]  = 13.28;
    args.params[CGCNeuron::delta_a] = 0.39;

    args.params[CGCNeuron::vh_b]    = -83.00;
    args.params[CGCNeuron::vs_b]    = -6.20;
    args.params[CGCNeuron::tau0_b]  = 266.75;
    args.params[CGCNeuron::delta_b] = 0.83;

    args.params[CGCNeuron::vh_n]    = -59.43;
    args.params[CGCNeuron::vs_n]    = 34.79;
    args.params[CGCNeuron::tau0_n]  = 14.52;
    args.params[CGCNeuron::delta_n] = 0.18;

    args.params[CGCNeuron::vh_e]    = -14.25;
    args.params[CGCNeuron::vs_e]    = 6.96;
    args.params[CGCNeuron::tau0_e]  = 3.81;
    args.params[CGCNeuron::delta_e] = 0.84;

    args.params[CGCNeuron::vh_f]    = -21.44;
    args.params[CGCNeuron::vs_f]    = -5.78;
    args.params[CGCNeuron::tau0_f]  = 34.68;
    args.params[CGCNeuron::delta_f] = 0.97;

    args.params[CGCNeuron::Vh_m] = -35.20;
    args.params[CGCNeuron::Vs_m] = 9.66;

    args.params[CGCNeuron::Vh_c] = -41.35;
    args.params[CGCNeuron::Vs_c] = 5.05;
    args.params[CGCNeuron::Vh_d] = -64.13;
    args.params[CGCNeuron::Vs_d] = -4.03;
}

// ─── Función principal: simular CGC con I_inj → devuelve spikes/min ──────────
/**
 * simulate_cgc(I_inj, verbose)
 *   Simula la CGC aislada durante SIM_TIME ms con corriente constante I_inj.
 *   Descarta los primeros TRANSIENT ms.
 *   Cuenta spikes por cruce ascendente del umbral V_THR.
 *   Retorna la tasa en spikes/min.
 *   Si verbose=true, imprime en stderr los detalles.
 */
static double simulate_cgc(double I_inj, bool verbose = false) {
    CGCNeuron::ConstructorArgs args;
    init_cgc(args);
    CGCNeuron cgc(args);

    // Condiciones iniciales (reposo ~-60 mV)
    cgc.set(CGCNeuron::v, -60.0);
    cgc.set(CGCNeuron::h, 1.0 / (1.0 + std::exp((-56.43 - (-60.0)) / -8.41)));
    cgc.set(CGCNeuron::r, 1.0 / (1.0 + std::exp((-47.03 - (-60.0)) / 20.55)));
    cgc.set(CGCNeuron::a, 1.0 / (1.0 + std::exp((-36.37 - (-60.0)) / 8.72)));
    cgc.set(CGCNeuron::b, 1.0 / (1.0 + std::exp((-83.00 - (-60.0)) / -6.20)));
    cgc.set(CGCNeuron::n, 1.0 / (1.0 + std::exp((-59.43 - (-60.0)) / 34.79)));
    cgc.set(CGCNeuron::e, 1.0 / (1.0 + std::exp((-14.25 - (-60.0)) / 6.96)));
    cgc.set(CGCNeuron::f, 1.0 / (1.0 + std::exp((-21.44 - (-60.0)) / -5.78)));

    int spike_count = 0;
    bool above_thr  = false;
    double v_prev   = -60.0;

    for (double t = 0.0; t < SIM_TIME; t += STEP) {
        cgc.add_synaptic_input(I_inj);
        cgc.step(STEP);

        double v = cgc.get(CGCNeuron::v);

        // Detección de spike: cruce ascendente del umbral
        if (t >= TRANSIENT) {
            if (!above_thr && v >= V_THR) {
                above_thr = true;
                ++spike_count;
            } else if (above_thr && v < V_THR) {
                above_thr = false;
            }
        }
        v_prev = v;
    }

    // Convertir a spikes/min (la ventana de conteo es SIM_WINDOW ms)
    double rate = spike_count * (60000.0 / SIM_WINDOW);  // spikes/min

    if (verbose) {
        std::cerr << "    I=" << std::fixed << std::setprecision(6) << I_inj
                  << " -> " << spike_count << " spikes en " << SIM_WINDOW/1000.0
                  << " s = " << std::setprecision(3) << rate << " spk/min" << std::endl;
    }

    return rate;
}

// ─── Búsqueda binaria ─────────────────────────────────────────────────────────
/**
 * binary_search_current(target_rate, I_lo, I_hi, verbose)
 *   Encuentra la corriente en [I_lo, I_hi] que produce target_rate spikes/min.
 *   Precondición: simulate_cgc(I_lo) <= target_rate <= simulate_cgc(I_hi)
 *                 (función monótona creciente con I).
 */
static double binary_search_current(double target_rate,
                                    double I_lo, double I_hi,
                                    bool verbose = false) {
    double rate_lo = simulate_cgc(I_lo, verbose);
    double rate_hi = simulate_cgc(I_hi, verbose);

    if (verbose) {
        std::cerr << "  Rango inicial: I_lo=" << I_lo << " (" << rate_lo
                  << " spk/min), I_hi=" << I_hi << " (" << rate_hi
                  << " spk/min), target=" << target_rate << std::endl;
    }

    // Verificar que el target está en el rango
    if (target_rate < rate_lo) {
        std::cerr << "[AVISO] La tasa objetivo (" << target_rate
                  << " spk/min) es MENOR que la tasa en I_lo=" << I_lo
                  << " (" << rate_lo << " spk/min). Devuelvo I_lo." << std::endl;
        return I_lo;
    }
    if (target_rate > rate_hi) {
        std::cerr << "[AVISO] La tasa objetivo (" << target_rate
                  << " spk/min) es MAYOR que la tasa en I_hi=" << I_hi
                  << " (" << rate_hi << " spk/min). Devuelvo I_hi." << std::endl;
        return I_hi;
    }

    double I_mid = I_lo;
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        I_mid = 0.5 * (I_lo + I_hi);

        if (std::fabs(I_hi - I_lo) < TOL_I) break;

        double rate_mid = simulate_cgc(I_mid, verbose);

        if (std::fabs(rate_mid - target_rate) < TOL_RATE) {
            if (verbose)
                std::cerr << "  Convergido en iteración " << iter+1 << std::endl;
            break;
        }

        // Ajustar intervalo: la frecuencia aumenta monótonamente con I
        if (rate_mid < target_rate)
            I_lo = I_mid;
        else
            I_hi = I_mid;
    }

    return I_mid;
}

// ─── main ─────────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    // Rangos de búsqueda (ajustables por línea de comandos)
    // I_min: corriente mínima (debería producir < 7 spk/min o silencio)
    // I_max: corriente máxima (debería producir > 60 spk/min)
    double I_search_min = 0.0;
    double I_search_max = 2.0;

    if (argc >= 3) {
        I_search_min = std::stod(argv[1]);
        I_search_max = std::stod(argv[2]);
    }

    bool verbose = false;
    if (argc >= 4 && std::string(argv[3]) == "--verbose")
        verbose = true;

    std::cerr << "============================================================" << std::endl;
    std::cerr << " cgc_freq_search: Búsqueda de corriente de estimulación CGC" << std::endl;
    std::cerr << "============================================================" << std::endl;
    std::cerr << " Ventana de simulación : " << SIM_TIME/1000.0 << " s"
              << " (transitorio=" << TRANSIENT/1000.0 << " s, conteo="
              << SIM_WINDOW/1000.0 << " s)" << std::endl;
    std::cerr << " Umbral de spike       : " << V_THR << " mV" << std::endl;
    std::cerr << " Paso de integración   : " << STEP << " ms" << std::endl;
    std::cerr << " Rango de búsqueda     : [" << I_search_min
              << ", " << I_search_max << "]" << std::endl;
    std::cerr << "------------------------------------------------------------" << std::endl;

    // ── Paso 1: Barrido inicial para acotar los rangos ──────────────────────
    // Comprobamos tasas en I_search_min e I_search_max para asegurarnos
    // de que los targets (7, 20, 42 y 60) están en el rango.
    std::cerr << "\n[1/5] Barrido inicial..." << std::endl;
    double rate_min = simulate_cgc(I_search_min, verbose);
    double rate_max = simulate_cgc(I_search_max, verbose);
    std::cerr << "  I=" << I_search_min << " -> " << rate_min << " spk/min" << std::endl;
    std::cerr << "  I=" << I_search_max << " -> " << rate_max << " spk/min" << std::endl;

    // ── Paso 2: Buscar corriente para 7 spk/min ──────────────────────────────
    std::cerr << "\n[2/5] Buscando corriente para 7 spk/min..." << std::endl;
    double I_7 = binary_search_current(7.0, I_search_min, I_search_max, verbose);
    double rate_7_found = simulate_cgc(I_7, false);
    std::cerr << "  >> I encontrada: " << std::fixed << std::setprecision(6) << I_7
              << " -> tasa verificada: " << std::setprecision(3) << rate_7_found
              << " spk/min" << std::endl;

    // ── Paso 3: Buscar corriente para 20 spk/min ─────────────────────────────
    std::cerr << "\n[3/5] Buscando corriente para 20 spk/min..." << std::endl;
    double I_20 = binary_search_current(20.0, I_search_min, I_search_max, verbose);
    double rate_20_found = simulate_cgc(I_20, false);
    std::cerr << "  >> I encontrada: " << std::fixed << std::setprecision(6) << I_20
              << " -> tasa verificada: " << std::setprecision(3) << rate_20_found
              << " spk/min" << std::endl;

    // ── Paso 4: Buscar corriente para 42 spk/min ─────────────────────────────
    std::cerr << "\n[4/5] Buscando corriente para 42 spk/min..." << std::endl;
    double I_42 = binary_search_current(42.0, I_search_min, I_search_max, verbose);
    double rate_42_found = simulate_cgc(I_42, false);
    std::cerr << "  >> I encontrada: " << std::fixed << std::setprecision(6) << I_42
              << " -> tasa verificada: " << std::setprecision(3) << rate_42_found
              << " spk/min" << std::endl;

    // ── Paso 5: Buscar corriente para 60 spk/min ─────────────────────────────
    std::cerr << "\n[5/5] Buscando corriente para 60 spk/min..." << std::endl;
    double I_60 = binary_search_current(60.0, I_search_min, I_search_max, verbose);
    double rate_60_found = simulate_cgc(I_60, false);
    std::cerr << "  >> I encontrada: " << std::fixed << std::setprecision(6) << I_60
              << " -> tasa verificada: " << std::setprecision(3) << rate_60_found
              << " spk/min" << std::endl;

    // ── Resultado final (a stdout, fácil de parsear) ─────────────────────────
    std::cerr << "\n============================================================" << std::endl;
    std::cerr << " RESULTADO FINAL" << std::endl;
    std::cerr << "============================================================" << std::endl;

    // stdout: formato limpio para scripts
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "# target_rate(spk/min)  I_corriente  tasa_verificada(spk/min)" << std::endl;
    std::cout << "7    " << I_7  << "  " << rate_7_found  << std::endl;
    std::cout << "20   " << I_20 << "  " << rate_20_found << std::endl;
    std::cout << "42   " << I_42 << "  " << rate_42_found << std::endl;
    std::cout << "60   " << I_60 << "  " << rate_60_found << std::endl;

    // Resumen legible también a stderr
    std::cerr << "  7  spk/min  ->  I = " << std::setprecision(6) << I_7
              << "  (verificado: " << std::setprecision(3) << rate_7_found  << " spk/min)" << std::endl;
    std::cerr << "  20 spk/min  ->  I = " << std::setprecision(6) << I_20
              << "  (verificado: " << std::setprecision(3) << rate_20_found << " spk/min)" << std::endl;
    std::cerr << "  42 spk/min  ->  I = " << std::setprecision(6) << I_42
              << "  (verificado: " << std::setprecision(3) << rate_42_found << " spk/min)" << std::endl;
    std::cerr << "  60 spk/min  ->  I = " << std::setprecision(6) << I_60
              << "  (verificado: " << std::setprecision(3) << rate_60_found << " spk/min)" << std::endl;

    return 0;
}
