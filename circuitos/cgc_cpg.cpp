/*************************************************************
 * cgc_cpg.cpp - CPG Completo + CGC (Cerebral Giant Cell)
 * 
 * Extiende el circuito CPG trifásico (SO, N1M, N2v, N3t)
 * añadiendo la neurona modulatoria CGC.
 * 
 * Conexiones CGC (Benjamin, 2012):
 *   CGC -----(excitación monosináptica)----> N1M
 *   CGC -----(excitación polisináptica)----> N2v
 *   CGC -----(excitación monosináptica)----> SO
 *   CGC -----(excitación monosináptica)----> N3t
 * 
 * La CGC es una neurona modulatoria que modula la actividad
 * del CPG de alimentación de Lymnaea stagnalis.
 * 
 *************************************************************/

#include <DifferentialNeuronWrapper.h>
#include <VavoulisModel.h>
#include <VavoulisCGCModel.h>
#include <GradualActivationSynapsis.h>
#include <RungeKutta4.h>
#include <SystemWrapper.h>
#include <iostream>

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<VavoulisModel<double>>, Integrator> Neuron;
typedef DifferentialNeuronWrapper<SystemWrapper<VavoulisCGCModel<double>>, Integrator> CGCNeuron;
typedef GradualActivationSynapsis<Neuron, Neuron, Integrator, double> Synapse;
typedef GradualActivationSynapsis<CGCNeuron, Neuron, Integrator, double> CGCSynapse;

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

  // CONFIGURACIÓN DE CGC (Cerebral Giant Cell - Neurona Modulatoria)
  // Modelo VavoulisCGCModel: 1 compartimento, 6 corrientes iónicas
  // Variables: v, h, r, a, b, n, e, f
  // Parámetros tomados de CGC.cpp (Tabla 1, probados)
  // Conexiones basadas en Benjamin (2012)
  CGCNeuron::ConstructorArgs args_cgc;

  // Capacitancia y potenciales de reversión (mV)
  args_cgc.params[CGCNeuron::cm]  = 1.0;      // Capacitancia de membrana (µF/cm²)
  args_cgc.params[CGCNeuron::vna] = 55.0;     // E_Na (mV)
  args_cgc.params[CGCNeuron::vk]  = -90.0;    // E_K (mV)
  args_cgc.params[CGCNeuron::vca] = 80.0;     // E_Ca (mV)

  // Conductancias máximas (mS/cm²)
  args_cgc.params[CGCNeuron::Gnat] = 1.68;    // Sodio transitorio
  args_cgc.params[CGCNeuron::Gnap] = 0.44;    // Sodio persistente
  args_cgc.params[CGCNeuron::Ga]   = 18.82;   // Potasio tipo A
  args_cgc.params[CGCNeuron::Gd]   = 1.20;    // Potasio rectificador retardado
  args_cgc.params[CGCNeuron::Glva] = 0.01;    // Calcio LVA (bajo umbral)
  args_cgc.params[CGCNeuron::Ghva] = 1.03;    // Calcio HVA (alto umbral)

  // h: inactivación Na transitorio
  args_cgc.params[CGCNeuron::vh_h]    = -56.43;  // V_half (mV)
  args_cgc.params[CGCNeuron::vs_h]    = -8.41;   // Pendiente (mV), negativa = inactivación
  args_cgc.params[CGCNeuron::tau0_h]  = 778.82;  // Constante de tiempo basal (ms)
  args_cgc.params[CGCNeuron::delta_h] = 0.03;    // Factor de asimetría de tau (adim.)

  // r: activación Na persistente
  args_cgc.params[CGCNeuron::vh_r]    = -47.03;  // V_half (mV)
  args_cgc.params[CGCNeuron::vs_r]    = 20.55;   // Pendiente (mV), positiva = activación
  args_cgc.params[CGCNeuron::tau0_r]  = 4.01;    // Constante de tiempo basal (ms)
  args_cgc.params[CGCNeuron::delta_r] = 1.00;    // Factor de asimetría de tau (adim.)

  // a: activación K tipo A
  args_cgc.params[CGCNeuron::vh_a]    = -36.37;  // V_half (mV)
  args_cgc.params[CGCNeuron::vs_a]    = 8.72;    // Pendiente (mV), positiva = activación
  args_cgc.params[CGCNeuron::tau0_a]  = 13.28;   // Constante de tiempo basal (ms)
  args_cgc.params[CGCNeuron::delta_a] = 0.39;    // Factor de asimetría de tau (adim.)

  // b: inactivación K tipo A
  args_cgc.params[CGCNeuron::vh_b]    = -83.00;  // V_half (mV)
  args_cgc.params[CGCNeuron::vs_b]    = -6.20;   // Pendiente (mV), negativa = inactivación
  args_cgc.params[CGCNeuron::tau0_b]  = 266.75;  // Constante de tiempo basal (ms)
  args_cgc.params[CGCNeuron::delta_b] = 0.83;    // Factor de asimetría de tau (adim.)

  // n: activación K rectificador retardado
  args_cgc.params[CGCNeuron::vh_n]    = -59.43;  // V_half (mV)
  args_cgc.params[CGCNeuron::vs_n]    = 34.79;   // Pendiente (mV), positiva = activación
  args_cgc.params[CGCNeuron::tau0_n]  = 14.52;   // Constante de tiempo basal (ms)
  args_cgc.params[CGCNeuron::delta_n] = 0.18;    // Factor de asimetría de tau (adim.)

  // e: activación Ca alto umbral (HVA)
  args_cgc.params[CGCNeuron::vh_e]    = -14.25;  // V_half (mV)
  args_cgc.params[CGCNeuron::vs_e]    = 6.96;    // Pendiente (mV), positiva = activación
  args_cgc.params[CGCNeuron::tau0_e]  = 3.81;    // Constante de tiempo basal (ms)
  args_cgc.params[CGCNeuron::delta_e] = 0.84;    // Factor de asimetría de tau (adim.)

  // f: inactivación Ca alto umbral (HVA)
  args_cgc.params[CGCNeuron::vh_f]    = -21.44;  // V_half (mV)
  args_cgc.params[CGCNeuron::vs_f]    = -5.78;   // Pendiente (mV), negativa = inactivación
  args_cgc.params[CGCNeuron::tau0_f]  = 34.68;   // Constante de tiempo basal (ms)
  args_cgc.params[CGCNeuron::delta_f] = 0.97;    // Factor de asimetría de tau (adim.)

  // m: activación Na transitorio (instantánea, sin dinámica propia)
  args_cgc.params[CGCNeuron::Vh_m] = -35.20;    // V_half (mV)
  args_cgc.params[CGCNeuron::Vs_m] = 9.66;      // Pendiente (mV)

  // c, d: activación/inactivación Ca bajo umbral LVA (instantáneas)
  args_cgc.params[CGCNeuron::Vh_c] = -41.35;    // V_half activación (mV)
  args_cgc.params[CGCNeuron::Vs_c] = 5.05;      // Pendiente activación (mV)
  args_cgc.params[CGCNeuron::Vh_d] = -64.13;    // V_half inactivación (mV)
  args_cgc.params[CGCNeuron::Vs_d] = -4.03;     // Pendiente inactivación (mV)

  CGCNeuron cgc(args_cgc);

  // Condiciones iniciales CGC (reposo basal aprox -60 mV)
  // x_inf(V) = 1/(1+exp((vh-V)/vs)), evaluadas en V = -60 mV
  cgc.set(CGCNeuron::v, -60.0);
  cgc.set(CGCNeuron::h, 1.0 / (1.0 + exp((-56.43 - (-60.0)) / -8.41)));
  cgc.set(CGCNeuron::r, 1.0 / (1.0 + exp((-47.03 - (-60.0)) / 20.55)));
  cgc.set(CGCNeuron::a, 1.0 / (1.0 + exp((-36.37 - (-60.0)) / 8.72)));
  cgc.set(CGCNeuron::b, 1.0 / (1.0 + exp((-83.00 - (-60.0)) / -6.20)));
  cgc.set(CGCNeuron::n, 1.0 / (1.0 + exp((-59.43 - (-60.0)) / 34.79)));
  cgc.set(CGCNeuron::e, 1.0 / (1.0 + exp((-14.25 - (-60.0)) / 6.96)));
  cgc.set(CGCNeuron::f, 1.0 / (1.0 + exp((-21.44 - (-60.0)) / -5.78)));

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

  // SINAPSIS DE CGC (Benjamin, 2012)
  // La CGC es una neurona modulatoria que usa sinapsis GradualActivation

  // Sinapsis CGC -> N1M (EXCITATORIA MONOSINÁPTICA)
  // CGC modula el CPG excitando N1M directamente (Benjamin, 2012)
  CGCSynapse::ConstructorArgs syn_cgc_to_n1m;
  syn_cgc_to_n1m.params[Synapse::esyn] = -9.0;        // Excitatoria 
  syn_cgc_to_n1m.params[Synapse::gsyn] = 3.0;        // Moderada para gating [cite: 474]
  syn_cgc_to_n1m.params[Synapse::tau_syn] = 175.0;   // Lenta (modulación tónica) 
  syn_cgc_to_n1m.params[Synapse::v_pre] = -67.0;
  syn_cgc_to_n1m.params[Synapse::v_r] = -40.0;
  syn_cgc_to_n1m.params[Synapse::dec_slope] = 2.5;
  CGCSynapse s_cgc_n1m(cgc, CGCNeuron::v, n1m, Neuron::v, syn_cgc_to_n1m, 1);

  // Sinapsis CGC -> N2v (EXCITATORIA POLISINÁPTICA)
  // Conexión indirecta/débil (Benjamin, 2012)
  CGCSynapse::ConstructorArgs syn_cgc_to_n2v;
  syn_cgc_to_n2v.params[Synapse::esyn] = -1.0;        // Excitatoria 
  syn_cgc_to_n2v.params[Synapse::gsyn] = 0.0;        // Necesaria para plateaus en N2v
  syn_cgc_to_n2v.params[Synapse::tau_syn] = 345.0;   
  syn_cgc_to_n2v.params[Synapse::v_pre] = -67.0;
  syn_cgc_to_n2v.params[Synapse::v_r] = -40.0;
  syn_cgc_to_n2v.params[Synapse::dec_slope] = 2.5;
  CGCSynapse s_cgc_n2v(cgc, CGCNeuron::v, n2v, Neuron::v, syn_cgc_to_n2v, 1);

  // Sinapsis CGC -> SO (EXCITATORIA MONOSINÁPTICA)
  // CGC excita al SO directamente (Benjamin, 2012)
  CGCSynapse::ConstructorArgs syn_cgc_to_so;
  syn_cgc_to_so.params[Synapse::esyn] = 10.0;         // Excitatoria 
  syn_cgc_to_so.params[Synapse::gsyn] = 0.0;         // Refuerza la capacidad del SO
  syn_cgc_to_so.params[Synapse::tau_syn] = 142.0;    
  syn_cgc_to_so.params[Synapse::v_pre] = -67.0;
  syn_cgc_to_so.params[Synapse::v_r] = -40.0;
  syn_cgc_to_so.params[Synapse::dec_slope] = 2.5;
  CGCSynapse s_cgc_so(cgc, CGCNeuron::v, so, Neuron::v, syn_cgc_to_so, 1);

  // Sinapsis CGC -> N3t (EXCITATORIA MONOSINÁPTICA)
  // CGC excita N3t directamente (Benjamin, 2012)
  CGCSynapse::ConstructorArgs syn_cgc_to_n3t;
  syn_cgc_to_n3t.params[Synapse::esyn] = -2.0;        // Excitatoria 
  syn_cgc_to_n3t.params[Synapse::gsyn] = 1.0;        
  syn_cgc_to_n3t.params[Synapse::tau_syn] = 164.0;   
  syn_cgc_to_n3t.params[Synapse::v_pre] = -67.0;
  syn_cgc_to_n3t.params[Synapse::v_r] = -40.0;
  syn_cgc_to_n3t.params[Synapse::dec_slope] = 2.5;
  CGCSynapse s_cgc_n3t(cgc, CGCNeuron::v, n3t, Neuron::v, syn_cgc_to_n3t, 1);

  // PARÁMETROS DE SIMULACIÓN
  const double step = 0.01;              // Paso de integración (ms)
  const double simulation_time = 50000;  // 100 segundos

  // Estimulación de SO para activar el CPG (como en Fig. 4C)
  const double t_stim_start = 100;       // Inicio del estímulo (ms)
  const double t_stim_end = 45000;        // Estímulo casi continuo
  
  const double I_drive_so = -8.5;       // Corriente despolarizante a SO
  const double I_drive_n1m = -6.0;      // Drive a N1M (más fuerte, es el "líder")
  const double I_drive_n2v = -2.0;      // Drive a N2v (más débil)
  const double I_drive_n3t = 0.0;       // Drive a N3t
  const double I_drive_cgc = 0.03;      // Drive a CGC 0.0298 es el minimo
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

    // Sinapsis de CGC (Benjamin, 2012)
    s_cgc_n1m.step(step);
    s_cgc_n2v.step(step);
    s_cgc_so.step(step);
    s_cgc_n3t.step(step);

    if (time >= t_stim_start && time <= t_stim_end) {
      so.add_synaptic_input(I_drive_so);
      n1m.add_synaptic_input(I_drive_n1m);  
      n2v.add_synaptic_input(I_drive_n2v);
      n3t.add_synaptic_input(I_drive_n3t);
      cgc.add_synaptic_input(I_drive_cgc);
    }

    // Integrar todas las neuronas
    n1m.step(step);
    n2v.step(step);
    n3t.step(step);
    so.step(step);
    cgc.step(step);

    // Salida: 29 columnas
    // tiempo(0), V_N1M_s(1), V_N1M_a(2), V_N2v_s(3), V_N2v_a(4), V_N3t_s(5), V_N3t_a(6), V_SO_s(7), V_SO_a(8), V_CGC(9),
    // I_n1m_n2v(10), I_n2v_n1m(11), I_n1m_n3t(12), I_n3t_n1m(13), I_n2v_n3t(14), I_n2v_so(15), I_so_n1m(16), I_so_n2v(17),
    // I_cgc_n1m(18), I_cgc_n2v(19), I_cgc_so(20), I_cgc_n3t(21),
    // p_N1M(22), p_N2v(23), q_N2v(24), p_N3t(25), q_N3t(26), p_SO(27), h_CGC(28)
    std::cout << time << " " 
              // Voltajes (9 valores: 8 del CPG + 1 de CGC)
              << n1m.get(Neuron::v) << " " 
              << n1m.get(Neuron::va) << " "
              << n2v.get(Neuron::v) << " " 
              << n2v.get(Neuron::va) << " "
              << n3t.get(Neuron::v) << " " 
              << n3t.get(Neuron::va) << " "
              << so.get(Neuron::v) << " " 
              << so.get(Neuron::va) << " "
              << cgc.get(CGCNeuron::v) << " "
              // Corrientes sinápticas CPG (8 valores)
              << s_n1m_n2v.get(Synapse::i) << " "
              << s_n2v_n1m.get(Synapse::i) << " "
              << s_n1m_n3t.get(Synapse::i) << " "
              << s_n3t_n1m.get(Synapse::i) << " "
              << s_n2v_n3t.get(Synapse::i) << " "
              << s_n2v_so.get(Synapse::i) << " "
              << s_so_n1m.get(Synapse::i) << " "
              << s_so_n2v.get(Synapse::i) << " "
              // Corrientes sinápticas CGC (4 valores)
              << s_cgc_n1m.get(CGCSynapse::i) << " "
              << s_cgc_n2v.get(CGCSynapse::i) << " "
              << s_cgc_so.get(CGCSynapse::i) << " "
              << s_cgc_n3t.get(CGCSynapse::i) << " "
              // Variables de gating (7 valores: 6 del CPG + 1 de CGC)
              << n1m.get(Neuron::p) << " "
              << n2v.get(Neuron::p) << " "
              << n2v.get(Neuron::q) << " "
              << n3t.get(Neuron::p) << " "
              << n3t.get(Neuron::q) << " "
              << so.get(Neuron::p) << " "
              << cgc.get(CGCNeuron::h)
              << std::endl;
  }

  return 0;
}