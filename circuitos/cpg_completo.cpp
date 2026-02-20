/*************************************************************
 * cpg_completo.cpp - CPG Completo de Alimentación de Lymnaea
 * 
 * Reproduce la Figura 4C del paper Vavoulis et al. (2007):
 * "Dynamic control of a central pattern generator circuit"
 * 
 * Red completa SO-driven:
 * 
 *        SO -----(excitación)----> N1M
 *        SO -----(excitación)----> N2v
 *        SO <----(inhibición)----- N2v
 * 
 *       N1M -----(excitación)----> N2v
 *       N1M <----(inhibición)----- N2v
 * 
 *       N1M -----(inhibición)----> N3t
 *       N1M <----(inhibición)----- N3t
 * 
 *       N2v -----(inhibición)----> N3t
 * 
 * El SO (Slow Oscillator) es una interneurona moduladora que:
 * - Controla la FRECUENCIA del ritmo de alimentación
 * - Es EXTRÍNSECA al CPG (no es necesaria para el patrón básico)
 * - Cuando está activa, produce ritmos de alta frecuencia (3-6s)
 * - Es inhibida durante la fase R por N2v
 * 
 * Parámetros sinápticos (Tabla 2, Vavoulis 2007):
 * 
 * SINAPSIS DEL CPG TRIFÁSICO:
 * - N1M -> N2v: Excitatoria lenta, g_syn = 0.077, E_syn = 0 mV, tau = 200 ms
 * - N2v -> N1M: Inhibitoria fuerte, g_syn = 50.0, E_syn = -90 mV, tau = 50 ms
 * - N1M -> N3t: Inhibitoria, g_syn = 5.0, E_syn = -90 mV, tau = 50 ms
 * - N3t -> N1M: Inhibitoria, g_syn = 8.0, E_syn = -90 mV, tau = 50 ms
 * - N2v -> N3t: Inhibitoria, g_syn = 2.0, E_syn = -90 mV, tau = 50 ms
 * 
 *************************************************************/

#include <DifferentialNeuronWrapper.h>
#include <VavoulisModel.h>
#include <GradualActivationSynapsis.h>
#include <RungeKutta4.h>
#include <SystemWrapper.h>
#include <iostream>

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<VavoulisModel<double>>, Integrator> Neuron;
typedef GradualActivationSynapsis<Neuron, Neuron, Integrator, double> Synapse;

int main(int argc, char **argv) {
  
  // CONFIGURACIÓN DE N1M (Interneurona de Protracción - Fase 1)
  Neuron::ConstructorArgs args_n1m;

  args_n1m.params[Neuron::n_type] = 1;     // n1m
  args_n1m.params[Neuron::tau_p] = 250.0;
  args_n1m.params[Neuron::tau_q] = 1.0;    // No se usa para N1M
  args_n1m.params[Neuron::g_eca] = 8.0;    // Acoplamiento axon->soma
  args_n1m.params[Neuron::g_ecs] = 8.0;    // Acoplamiento soma->axon

  Neuron n1m(args_n1m);
  
  n1m.set(Neuron::v, -67);
  n1m.set(Neuron::va, -67);
  n1m.set(Neuron::p, 1 / (1 + exp((-38.8 - (-67.0)) / 10.0)));
  n1m.set(Neuron::q, 0.0);
  n1m.set(Neuron::h, 1 / (1 + exp((-55.2 - (-67.0))/-7.1)));
  n1m.set(Neuron::n, 1.0 / (1.0 + exp((-30.0 - (-67.0)) / 17.4)));

  // CONFIGURACIÓN DE N2v (Interneurona de Rasp - Fase 2)
  Neuron::ConstructorArgs args_n2v;

  args_n2v.params[Neuron::n_type] = 2;     // n2v
  args_n2v.params[Neuron::tau_p] = 1.0;    // No usado para N2v
  args_n2v.params[Neuron::tau_q] = 1.0;    // No usado para N2v
  args_n2v.params[Neuron::g_eca] = 0.06;   // Acoplamiento axon->soma
  args_n2v.params[Neuron::g_ecs] = 0.55;   // Acoplamiento soma->axon

  Neuron n2v(args_n2v);

  // Condiciones iniciales (estado de reposo)
  n2v.set(Neuron::v, -67);
  n2v.set(Neuron::va, -67);
  n2v.set(Neuron::p, 1 / (1 + exp((-51 - (-67.0))/10.3)));
  n2v.set(Neuron::q, 1 / (1 + exp((-45 - (-67.0))/-3)));
  n2v.set(Neuron::h, 1 / (1 + exp((-55.2 - (-67.0))/-7.1)));
  n2v.set(Neuron::n, 1.0 / (1.0 + exp((-30.0 - (-67.0)) / 17.4)));

  // CONFIGURACIÓN DE N3t (Interneurona de Swallow - Fase 3)
  Neuron::ConstructorArgs args_n3t;

  args_n3t.params[Neuron::n_type]  = 3;     // n3t
  args_n3t.params[Neuron::tau_p] = 4.0;
  args_n3t.params[Neuron::tau_q] = 400.0;
  args_n3t.params[Neuron::g_eca] = 8.0;    // Acoplamiento axon->soma
  args_n3t.params[Neuron::g_ecs] = 8.0;    // Acoplamiento soma->axon

  Neuron n3t(args_n3t);
  
  n3t.set(Neuron::v, -67.0);
  n3t.set(Neuron::va, -67.0);
  n3t.set(Neuron::p, 1 / (1 + exp((-61.6 - (-67))/5.6)));   
  n3t.set(Neuron::q, 1 / (1 + exp((-73.2 - (-67))/-5.1)));   
  n3t.set(Neuron::h, 1 / (1 + exp((-55.2 - (-67.0))/-7.1)));
  n3t.set(Neuron::n, 1.0 / (1.0 + exp((-30.0 - (-67.0)) / 17.4)));


  // CONFIGURACIÓN DE SO (Slow Oscillator - Modulación de Frecuencia)
  Neuron::ConstructorArgs args_so;

  args_so.params[Neuron::n_type] = 0;
  args_so.params[Neuron::tau_p] = 1.0;    // No usado
  args_so.params[Neuron::tau_q] = 1.0;    // No usado
  args_so.params[Neuron::g_eca] = 8.0;    // Acoplamiento axon->soma
  args_so.params[Neuron::g_ecs] = 8.0;    // Acoplamiento soma->axon

  Neuron so(args_so);
  
  // Condiciones iniciales para SO
  so.set(Neuron::v, -67);
  so.set(Neuron::va, -67);
  so.set(Neuron::p, 0.0);    // No usado para SO
  so.set(Neuron::q, 0.0);    // No usado para SO
  so.set(Neuron::h, 1 / (1 + exp((-55.2 - (-67.0))/-7.1)));
  so.set(Neuron::n, 1.0 / (1.0 + exp((-30.0 - (-67.0)) / 17.4)));

  // SINAPSIS DEL CPG TRIFÁSICO
  
  // Sinapsis N1M -> N2v (EXCITATORIA LENTA)
  Synapse::ConstructorArgs syn_n1m_to_n2v;
  syn_n1m_to_n2v.params[Synapse::esyn] = 0.0;        // E_syn = 0 mV (excitatoria)
  syn_n1m_to_n2v.params[Synapse::gsyn] = 0.077;      // Conductancia (Tabla 2)
  syn_n1m_to_n2v.params[Synapse::tau_syn] = 200.0;   // tau = 200 ms (lenta)
  syn_n1m_to_n2v.params[Synapse::v_pre] = -67.0;
  syn_n1m_to_n2v.params[Synapse::v_r] = -40.0;
  syn_n1m_to_n2v.params[Synapse::dec_slope] = 2.5;
  Synapse s_n1m_n2v(n1m, Neuron::v, n2v, Neuron::v, syn_n1m_to_n2v, 1);

  // Sinapsis N2v -> N1M (INHIBITORIA FUERTE)
  Synapse::ConstructorArgs syn_n2v_to_n1m;
  syn_n2v_to_n1m.params[Synapse::esyn] = -90.0;      // E_syn = -90 mV
  syn_n2v_to_n1m.params[Synapse::gsyn] = 50.0;
  syn_n2v_to_n1m.params[Synapse::tau_syn] = 50.0;    // tau = 50 ms
  syn_n2v_to_n1m.params[Synapse::v_pre] = -67.0;
  syn_n2v_to_n1m.params[Synapse::v_r] = -40.0;
  syn_n2v_to_n1m.params[Synapse::dec_slope] = 2.5;
  Synapse s_n2v_n1m(n2v, Neuron::v, n1m, Neuron::v, syn_n2v_to_n1m, 1);

  // Sinapsis N1M -> N3t (INHIBITORIA)
  Synapse::ConstructorArgs syn_n1m_to_n3t;
  syn_n1m_to_n3t.params[Synapse::esyn] = -90.0;
  syn_n1m_to_n3t.params[Synapse::gsyn] = 0.5;        // g_syn = 0.5 (Tabla 2)
  syn_n1m_to_n3t.params[Synapse::tau_syn] = 50.0;
  syn_n1m_to_n3t.params[Synapse::v_pre] = -67.0;
  syn_n1m_to_n3t.params[Synapse::v_r] = -40.0;
  syn_n1m_to_n3t.params[Synapse::dec_slope] = 2.5;
  Synapse s_n1m_n3t(n1m, Neuron::v, n3t, Neuron::v, syn_n1m_to_n3t, 1);

  // Sinapsis N3t -> N1M (INHIBITORIA)
  Synapse::ConstructorArgs syn_n3t_to_n1m;
  syn_n3t_to_n1m.params[Synapse::esyn] = -90.0;
  syn_n3t_to_n1m.params[Synapse::gsyn] = 8.0;        // g_syn = 8.0 (Tabla 2)
  syn_n3t_to_n1m.params[Synapse::tau_syn] = 50.0;
  syn_n3t_to_n1m.params[Synapse::v_pre] = -67.0;
  syn_n3t_to_n1m.params[Synapse::v_r] = -40.0;
  syn_n3t_to_n1m.params[Synapse::dec_slope] = 2.5;
  Synapse s_n3t_n1m(n3t, Neuron::v, n1m, Neuron::v, syn_n3t_to_n1m, 1);

  // Sinapsis N2v -> N3t (INHIBITORIA)
  Synapse::ConstructorArgs syn_n2v_to_n3t;
  syn_n2v_to_n3t.params[Synapse::esyn] = -90.0;
  syn_n2v_to_n3t.params[Synapse::gsyn] = 2.0;        // g_syn = 2.0 (Tabla 2)
  syn_n2v_to_n3t.params[Synapse::tau_syn] = 50.0;
  syn_n2v_to_n3t.params[Synapse::v_pre] = -67.0;
  syn_n2v_to_n3t.params[Synapse::v_r] = -40.0;
  syn_n2v_to_n3t.params[Synapse::dec_slope] = 2.5;
  Synapse s_n2v_n3t(n2v, Neuron::v, n3t, Neuron::v, syn_n2v_to_n3t, 1);

  // SINAPSIS CON SO (Tabla 2, Vavoulis 2007)

  // Sinapsis N2v -> SO (INHIBITORIA) - Inhibe SO durante fase R
  Synapse::ConstructorArgs syn_n2v_to_so;
  syn_n2v_to_so.params[Synapse::esyn] = -90.0;       // E_syn = -90 mV (inhibitoria)
  syn_n2v_to_so.params[Synapse::gsyn] = 8.0;         // g_syn = 8.0
  syn_n2v_to_so.params[Synapse::tau_syn] = 50.0;     // tau = 50 ms
  syn_n2v_to_so.params[Synapse::v_pre] = -67.0;
  syn_n2v_to_so.params[Synapse::v_r] = -40.0;
  syn_n2v_to_so.params[Synapse::dec_slope] = 2.5;
  Synapse s_n2v_so(n2v, Neuron::v, so, Neuron::v, syn_n2v_to_so, 1);

  // Sinapsis SO -> N1M (EXCITATORIA LENTA) - Activa el CPG
  Synapse::ConstructorArgs syn_so_to_n1m;
  syn_so_to_n1m.params[Synapse::esyn] = 0.0;         // E_syn = 0 mV (excitatoria)
  syn_so_to_n1m.params[Synapse::gsyn] = 4.0;         // g_syn = 4.0
  syn_so_to_n1m.params[Synapse::tau_syn] = 200.0;    // tau = 200 ms (lenta)
  syn_so_to_n1m.params[Synapse::v_pre] = -67.0;
  syn_so_to_n1m.params[Synapse::v_r] = -40.0;
  syn_so_to_n1m.params[Synapse::dec_slope] = 2.5;
  Synapse s_so_n1m(so, Neuron::v, n1m, Neuron::v, syn_so_to_n1m, 1);

  // Sinapsis SO -> N2v (EXCITATORIA LENTA) - Acelera activación de N2v
  Synapse::ConstructorArgs syn_so_to_n2v;
  syn_so_to_n2v.params[Synapse::esyn] = 0.0;         // E_syn = 0 mV (excitatoria)
  syn_so_to_n2v.params[Synapse::gsyn] = 1.0;         // g_syn = 1.0
  syn_so_to_n2v.params[Synapse::tau_syn] = 200.0;    // tau = 200 ms (lenta)
  syn_so_to_n2v.params[Synapse::v_pre] = -67.0;
  syn_so_to_n2v.params[Synapse::v_r] = -40.0;
  syn_so_to_n2v.params[Synapse::dec_slope] = 2.5;
  Synapse s_so_n2v(so, Neuron::v, n2v, Neuron::v, syn_so_to_n2v, 1);

  // PARÁMETROS DE SIMULACIÓN
  const double step = 0.01;              // Paso de integración (ms)
  const double simulation_time = 10000;  // 10 segundos

  // Estimulación de SO para activar el CPG (como en Fig. 4C)
  const double t_stim_start = 100;       // Inicio del estímulo (ms)
  const double t_stim_end = 9500;        // Estímulo casi continuo
  
  // Drive a SO: controla la frecuencia del ritmo
  // Valores mayores -> ritmos más rápidos (hasta ~0.33 Hz)
  const double I_drive_so = -8.5;       // Corriente despolarizante a SO
  const double I_drive_n1m = -6.0;       // Drive a N1M (más fuerte, es el "líder")
  const double I_drive_n2v = -2.0;       // Drive a N2v (más débil)
  const double I_drive_n3t = 0.0;       // Drive a N3t
  // BUCLE DE SIMULACIÓN
  for (double time = 0; time < simulation_time; time += step) {
    
    // Actualizar todas las sinapsis (8 en total)
    // CPG trifásico
    s_n1m_n2v.step(step);
    s_n2v_n1m.step(step);
    s_n1m_n3t.step(step);
    s_n3t_n1m.step(step);
    s_n2v_n3t.step(step);

    // Sinapsis con SO
    s_n2v_so.step(step);
    s_so_n1m.step(step);
    s_so_n2v.step(step);

    if (time >= t_stim_start && time <= t_stim_end) {
      so.add_synaptic_input(I_drive_so);
      n1m.add_synaptic_input(I_drive_n1m);  
      n2v.add_synaptic_input(I_drive_n2v);
      n3t.add_synaptic_input(I_drive_n3t);
    }

    // Integrar todas las neuronas
    n1m.step(step);
    n2v.step(step);
    n3t.step(step);
    so.step(step);

    // Salida: 23 columnas
    // tiempo, V_N1M_s, V_N1M_a, V_N2v_s, V_N2v_a, V_N3t_s, V_N3t_a, V_SO_s, V_SO_a,
    // I_n1m_n2v, I_n2v_n1m, I_n1m_n3t, I_n3t_n1m, I_n2v_n3t, I_n2v_so, I_so_n1m, I_so_n2v,
    // p_N1M, p_N2v, q_N2v, p_N3t, q_N3t, p_SO
    std::cout << time << " " 
              // Voltajes (8 valores)
              << n1m.get(Neuron::v) << " " 
              << n1m.get(Neuron::va) << " "
              << n2v.get(Neuron::v) << " " 
              << n2v.get(Neuron::va) << " "
              << n3t.get(Neuron::v) << " " 
              << n3t.get(Neuron::va) << " "
              << so.get(Neuron::v) << " " 
              << so.get(Neuron::va) << " "
              // Corrientes sinápticas (8 valores)
              << s_n1m_n2v.get(Synapse::i) << " "
              << s_n2v_n1m.get(Synapse::i) << " "
              << s_n1m_n3t.get(Synapse::i) << " "
              << s_n3t_n1m.get(Synapse::i) << " "
              << s_n2v_n3t.get(Synapse::i) << " "
              << s_n2v_so.get(Synapse::i) << " "
              << s_so_n1m.get(Synapse::i) << " "
              << s_so_n2v.get(Synapse::i) << " "
              // Variables de gating (6 valores)
              << n1m.get(Neuron::p) << " "
              << n2v.get(Neuron::p) << " "
              << n2v.get(Neuron::q) << " "
              << n3t.get(Neuron::p) << " "
              << n3t.get(Neuron::q) << " "
              << so.get(Neuron::p)
              << std::endl;
  }

  return 0;
}
