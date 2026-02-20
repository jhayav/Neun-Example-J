/*************************************************************
 * n1-2-3.cpp - Circuito N1M <--> N2v (Fase trifásica)
 *                            \/
 *                           N3t
 * Reproduce la Figura 4B del paper Vavoulis et al. (2007):
 * "Dynamic control of a central pattern generator circuit"
 * 
 * El circuito N1M <--> N2v 
 *                  \/
 *                  N3t
 * forma un CPG de alimentación trifásico en el que:
 * 
 *   N1M ----(excitación)----> N2v
 *   N1M <---(inhibición)----- N2v
 * 
 *   N1M ----(inhibición)----> N3t
 *   N1M <---(inhibición)----- N3t
 * 
 *   N2v ----(inhibición)----> N3t
 *   
 * 
 * Parámetros sinápticos (Tabla 2, Vavoulis 2007):
 * - N1M -> N2v: Excitatoria lenta, g_syn = 0.077, E_syn = 0 mV, tau = 200 ms
 * - N2v -> N1M: Inhibitoria fuerte, g_syn = 50.0, E_syn = -90 mV, tau = 50 ms
 * - N1M -> N3t: Inhibitoria moderada, g_syn = 5.0, E_syn = -90 mV, tau = 50 ms
 * - N3t -> N1M: Inhibitoria débil, g_syn = 8.0, E_syn = -90 mV, tau = 50 ms
 * - N2v -> N3t: Inhibitoria moderada, g_syn = 2.0, E_syn = -90 mV, tau = 50 ms
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
  n2v.set(Neuron::v, -67);    // Reposo cerca de E_leak
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


  // SINAPSIS

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

  const double step = 0.01;              // Paso de integración (ms)
  const double simulation_time = 10000;  // 10 segundos para ver oscilaciones

  // Drive tónico continuo a ambas neuronas
  // Simula la entrada excitatoria de SO que activa el sistema
  const double t_stim_start = 100;       // Inicio del estímulo (ms)
  const double t_stim_end = 9500;
  
  const double I_drive_n1m = -6.0;
  const double I_drive_n2v = -1.0;
  const double I_drive_n3t = -3.0;

  for (double time = 0; time < simulation_time; time += step) {
    
    // Actualizar sinapsis con activación gradual - LAS 5 sinapsis de la Tabla 2
    s_n1m_n2v.step(step);
    s_n2v_n1m.step(step);
    s_n1m_n3t.step(step);
    s_n3t_n1m.step(step);
    s_n2v_n3t.step(step);

    // Drive tónico diferencial a cada neurona (simula entrada de SO)
    // Este drive es necesario para que el half-center oscillator funcione
    if (time >= t_stim_start && time <= t_stim_end) {
      n1m.add_synaptic_input(I_drive_n1m);
      n2v.add_synaptic_input(I_drive_n2v);
      n3t.add_synaptic_input(I_drive_n3t);
    }

    n1m.step(step);
    n2v.step(step);
    n3t.step(step);

    // Salida: tiempo, V_N1M_soma, V_N1M_axon, V_N2v_soma, V_N2v_axon, V_N3t_soma, V_N3t_axon,
    //         I_syn_N1M->N2v, I_syn_N2v->N1M, I_syn_N1M->N3t, I_syn_N3t->N1M, I_syn_N2v->N3t,
    //         p_N1M, p_N2v, q_N2v, p_N3t, q_N3t
    std::cout << time << " " 
              << n1m.get(Neuron::v) << " " 
              << n1m.get(Neuron::va) << " "
              << n2v.get(Neuron::v) << " " 
              << n2v.get(Neuron::va) << " "
              << n3t.get(Neuron::v) << " " 
              << n3t.get(Neuron::va) << " "
              << s_n1m_n2v.get(Synapse::i) << " "
              << s_n2v_n1m.get(Synapse::i) << " "
              << s_n1m_n3t.get(Synapse::i) << " "
              << s_n3t_n1m.get(Synapse::i) << " "
              << s_n2v_n3t.get(Synapse::i) << " "
              << n1m.get(Neuron::p) << " "
              << n2v.get(Neuron::p) << " "
              << n2v.get(Neuron::q) << " "
              << n3t.get(Neuron::p) << " "
              << n3t.get(Neuron::q)
              << std::endl;
  }

  return 0;
}
