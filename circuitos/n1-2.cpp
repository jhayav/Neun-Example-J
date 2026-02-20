/*************************************************************
 * n1-2.cpp - Circuito N1M <-> N2v (Core Pacemaker)
 * 
 * Reproduce la Figura 4A(ii) del paper Vavoulis et al. (2007):
 * "Dynamic control of a central pattern generator circuit"
 * 
 * El circuito N1M <-> N2v forma el "core pacemaker" del CPG de
 * alimentación de Lymnaea stagnalis. Las conexiones son:
 * 
 *   N1M ----(excitación)----> N2v
 *   N1M <---(inhibición)----- N2v
 * 
 * Ambas sinapsis son QUÍMICAS (no eléctricas) según Tabla 2 del paper.
 * 
 * IMPORTANTE: El paper usa ecuaciones diferenciales de SEGUNDO ORDEN
 * para las sinapsis (Ec. 3 y 4, Vavoulis 2007):
 * 
 *   tau_syn * dr/dt = r_inf - r
 *   tau_syn * ds/dt = r - s
 *   I_syn = g_syn * s * (V_post - E_syn)
 * 
 * Donde r_inf = 1 / (1 + exp((V_r - V_pre) / dec_slope))
 * 
 * Por ello usamos GradualActivationSynapsis que implementa este modelo.
 * 
 * Parámetros sinápticos (Tabla 2, Vavoulis 2007):
 * - N1M -> N2v: Excitatoria lenta, g_syn = 0.077, E_syn = 0 mV, tau = 200 ms
 * - N2v -> N1M: Inhibitoria fuerte, g_syn = 50.0, E_syn = -90 mV, tau = 50 ms
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
  
  // CONFIGURACIÓN DE N1M (Interneurona de Protracción)
  Neuron::ConstructorArgs args_n1m;

  args_n1m.params[Neuron::n_type] = 1;     // n1m
  args_n1m.params[Neuron::tau_p] = 250.0;  // Constante de tiempo lenta (ms)
  args_n1m.params[Neuron::tau_q] = 1.0;    // No se usa para N1M
  args_n1m.params[Neuron::g_eca] = 8.0;    // Acoplamiento axon->soma
  args_n1m.params[Neuron::g_ecs] = 8.0;    // Acoplamiento soma->axon

  Neuron n1m(args_n1m);
  
  // Condiciones iniciales (estado de reposo)
  n1m.set(Neuron::v, -67);
  n1m.set(Neuron::va, -67);
  n1m.set(Neuron::p, 1 / (1 + exp((-38.8 - (-67.0)) / 10.0)));
  n1m.set(Neuron::q, 0.0);
  n1m.set(Neuron::h, 1 / (1 + exp((-55.2 - (-67.0))/-7.1)));
  n1m.set(Neuron::n, 1.0 / (1.0 + exp((-30.0 - (-67.0)) / 17.4)));

  // CONFIGURACIÓN DE N2v (Interneurona de Rasp)
  Neuron::ConstructorArgs args_n2v;

  args_n2v.params[Neuron::n_type] = 2;     // n2v
  args_n2v.params[Neuron::tau_p] = 1.0;    // No usado para N2v
  args_n2v.params[Neuron::tau_q] = 1.0;    // No usado para N2v
  args_n2v.params[Neuron::g_eca] = 0.06;   // Acoplamiento axon->soma (muy debil)
  args_n2v.params[Neuron::g_ecs] = 0.55;   // Acoplamiento soma->axon
  
  Neuron n2v(args_n2v);

  // Condiciones iniciales (estado de reposo)
  n2v.set(Neuron::v, -67);
  n2v.set(Neuron::va, -67);
  n2v.set(Neuron::p, 1 / (1 + exp((-51 - (-67.0))/10.3)));
  n2v.set(Neuron::q, 1 / (1 + exp((-45 - (-67.0))/-3)));
  n2v.set(Neuron::h, 1 / (1 + exp((-55.2 - (-67.0))/-7.1)));
  n2v.set(Neuron::n, 1.0 / (1.0 + exp((-30.0 - (-67.0)) / 17.4)));

  // SINAPSIS QUÍMICAS CON ACTIVACIÓN GRADUAL (Tabla 2, Vavoulis 2007)
  // Modelo de segundo orden: tau*dr/dt = r_inf - r, tau*ds/dt = r - s
  // 
  // GradualActivationSynapsisModel parámetros:
  //   enum parameter { esyn, gsyn, tau_syn, v_pre, v_r, dec_slope, n_parameters };
  
  // Sinapsis N1M -> N2v (EXCITATORIA LENTA)
  // Cuando N1M dispara, excita gradualmente a N2v para iniciar la fase de rasp
  Synapse::ConstructorArgs syn_n1m_to_n2v;
  
  syn_n1m_to_n2v.params[Synapse::esyn] = 0.0;        // E_syn = 0 mV (excitatoria)
  syn_n1m_to_n2v.params[Synapse::gsyn] = 0.077;      // Conductancia
  syn_n1m_to_n2v.params[Synapse::tau_syn] = 200.0;   // tau = 200 ms
  syn_n1m_to_n2v.params[Synapse::v_pre] = -67.0;     // Se actualiza dinámicamente
  syn_n1m_to_n2v.params[Synapse::v_r] = -40.0;       // Umbral de activación
  syn_n1m_to_n2v.params[Synapse::dec_slope] = 2.5;   // Pendiente de la sigmoide

  Synapse s_n1m_n2v(n1m, Neuron::v, n2v, Neuron::v, syn_n1m_to_n2v, 1);

  // Sinapsis N2v -> N1M (INHIBITORIA FUERTE)
  // Cuando N2v dispara, inhibe fuertemente a N1M para terminar la protracción
  Synapse::ConstructorArgs syn_n2v_to_n1m;
  
  syn_n2v_to_n1m.params[Synapse::esyn] = -90.0;      // E_syn = -90 mV (inhibitoria)
  syn_n2v_to_n1m.params[Synapse::gsyn] = 50.0;       // Conductancia
  syn_n2v_to_n1m.params[Synapse::tau_syn] = 50.0;    // tau = 50 ms
  syn_n2v_to_n1m.params[Synapse::v_pre] = -67.0;     // Se actualiza dinámicamente
  syn_n2v_to_n1m.params[Synapse::v_r] = -40.0;       // Umbral de activación
  syn_n2v_to_n1m.params[Synapse::dec_slope] = 2.5;   // Pendiente de la sigmoide

  Synapse s_n2v_n1m(n2v, Neuron::v, n1m, Neuron::v, syn_n2v_to_n1m, 1);

  const double step = 0.01;              // Paso de integración (ms)
  const double simulation_time = 10000;

  // Drive tónico continuo a ambas neuronas
  // Simula la entrada excitatoria de SO que activa el sistema
  const double t_stim_start = 100;       // Inicio del estímulo (ms)
  const double t_stim_end = 9500;        // Estímulo casi continuo
  const double I_drive_n1m = -6.0;       // Drive a N1M
  const double I_drive_n2v = -1.5;       // Drive a N2v

  for (double time = 0; time < simulation_time; time += step) {
    
    // Actualizar sinapsis con activación gradual
    // NOTA: step(h) actualiza v_pre internamente y aplica I_syn a la neurona postsináptica
    s_n1m_n2v.step(step);  // N1M -> N2v
    s_n2v_n1m.step(step);  // N2v -> N1M

    if (time >= t_stim_start && time <= t_stim_end) {
      n1m.add_synaptic_input(I_drive_n1m);  // Drive más fuerte a N1M
      n2v.add_synaptic_input(I_drive_n2v);  // Drive más débil a N2v
    }

    n1m.step(step);
    n2v.step(step);

    // Salida: tiempo, V_N1M_soma, V_N1M_axon, V_N2v_soma, V_N2v_axon, 
    //         I_syn_N1M->N2v, I_syn_N2v->N1M, p_N1M, p_N2v, q_N2v
    std::cout << time << " " 
              << n1m.get(Neuron::v) << " " 
              << n1m.get(Neuron::va) << " "
              << n2v.get(Neuron::v) << " " 
              << n2v.get(Neuron::va) << " "
              << s_n1m_n2v.get(Synapse::i) << " "
              << s_n2v_n1m.get(Synapse::i) << " "
              << n1m.get(Neuron::p) << " "
              << n2v.get(Neuron::p) << " "
              << n2v.get(Neuron::q)
              << std::endl;
  }

  return 0;
}
