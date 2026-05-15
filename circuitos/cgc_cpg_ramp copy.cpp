/*************************************************************
 * cgc_cpg_ramp.cpp - CPG+CGC con protocolo de rampa de corriente
 * 
 * Uso: ./cgc_ramp <NEURON>
 *   NEURON = N1M | N3t | SO | CGC
 *
 *************************************************************/

#include <DifferentialNeuronWrapper.h>
#include <VavoulisModel.h>
#include <VavoulisCGCModel.h>
#include <GradualActivationSynapsis.h>
#include <RungeKutta4.h>
#include <SystemWrapper.h>
#include <iostream>
#include <cmath>
#include <cstring>

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<VavoulisModel<double>>, Integrator> Neuron;
typedef DifferentialNeuronWrapper<SystemWrapper<VavoulisCGCModel<double>>, Integrator> CGCNeuron;
typedef GradualActivationSynapsis<Neuron, Neuron, Integrator, double> Synapse;
typedef GradualActivationSynapsis<CGCNeuron, Neuron, Integrator, double> CGCSynapse;

// Rampa triangular generalizada.
// Interpola linealmente entre c_min y c_max en n_steps pasos.
// 1 triángulo completo: c_min->c_max->c_min
double compute_ramp(double time, double t_start,
                    double c_min, double c_max, int n_steps,
                    double interval) {
    if (time < t_start) return c_min;
    double t_rel = time - t_start;
    double half = n_steps * interval;
    double full = 2.0 * half;           // 1 triángulo = subida + bajada
    if (t_rel >= full) return c_min;
    int idx;
    if (t_rel < half) {
        idx = (int)(t_rel / interval);
        if (idx > n_steps) idx = n_steps;
    } else {
        idx = n_steps - (int)((t_rel - half) / interval);
        if (idx < 0) idx = 0;
    }
    return c_min + (c_max - c_min) * (double)idx / (double)n_steps;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Uso: " << argv[0] << " <N1M|N3t|SO|CGC>" << std::endl;
    return 1;
  }
  const char* target = argv[1];

  // --- Parámetros según Table 1 (con signo directo como en cgc_cpg.cpp) ---
  // Negativo = despolarizante para neuronas CPG (SO, N1M, N2v, N3t)
  // Positivo = despolarizante para CGC
  double I_so, I_n1m, I_n2v, I_n3t, I_cgc;
  double c_min, c_max;
  int n_ramp_steps;
  const double ramp_interval = 4600.0;
  int ramp_target = -1; // 0=SO, 1=N1M, 3=N3t, 4=CGC

if (strcmp(target, "N1M") == 0) {
    // Table 1: SO=-8.5, N1M=rampa, N2v=-2, N3t=0 | rampa: 0→-10.5 (21 steps)
    I_so = -8.5;  I_n1m = 0.0;   I_n2v = -2.0;  I_n3t = 0.0;
    I_cgc = 0.03;
    c_min = 0.0;  c_max = -10.5; n_ramp_steps = 21;
    ramp_target = 1;
  } else if (strcmp(target, "N3t") == 0) {
    // Table 1 adaptada: SO=-9, N1M=-10, N2v=-1, N3t=rampa | rampa: 0→-5 (20 steps)
    I_so = -9.0;  I_n1m = -10.0; I_n2v = -1.0;  I_n3t = 0.0;
    I_cgc = 0.03;
    c_min = 0.0;  c_max = -5.0;  n_ramp_steps = 20;
    ramp_target = 3;
  } else if (strcmp(target, "SO") == 0) {
    // Table 1 adaptada: SO=rampa, N1M=-10, N2v=-1, N3t=-4 | rampa: -8.2 -13 (20 steps) 
    I_so = 0.0;   I_n1m = -10.0; I_n2v = -1.0;  I_n3t = -4.0;
    I_cgc = 0.03;
    c_min = -8.2; c_max = -13.0; n_ramp_steps = 20; 
    ramp_target = 0;
  } else if (strcmp(target, "CGC") == 0) {
    // Nuevo: base del circuito cgc_cpg.cpp, rampa en CGC: 0→0.15 (20 steps)
    I_so = -8.5;  I_n1m = -6.0;  I_n2v = -2.0;  I_n3t = 0.0;
    I_cgc = 0.0;
    c_min = 0.032959;  c_max = 0.05;  n_ramp_steps = 20;
    // c_max = 0.042969;
    ramp_target = 4;
  } else {
    std::cerr << "Neurona no válida: " << target << std::endl;
    return 1;
  }

  // =================== NEURONAS ===================
  Neuron::ConstructorArgs args_n1m;
  args_n1m.params[Neuron::n_type] = 1;
  args_n1m.params[Neuron::tau_p] = 250.0;
  args_n1m.params[Neuron::tau_q] = 1.0;
  args_n1m.params[Neuron::g_eca] = 8.0;
  args_n1m.params[Neuron::g_ecs] = 8.0;
  Neuron n1m(args_n1m);
  n1m.set(Neuron::v, -67); n1m.set(Neuron::va, -67);
  n1m.set(Neuron::p, 1/(1+exp((-38.8-(-67.0))/10.0)));
  n1m.set(Neuron::q, 0.0);
  n1m.set(Neuron::h, 1/(1+exp((-55.2-(-67.0))/-7.1)));
  n1m.set(Neuron::n, 1.0/(1.0+exp((-30.0-(-67.0))/17.4)));

  Neuron::ConstructorArgs args_n2v;
  args_n2v.params[Neuron::n_type] = 2;
  args_n2v.params[Neuron::tau_p] = 1.0;
  args_n2v.params[Neuron::tau_q] = 1.0;
  args_n2v.params[Neuron::g_eca] = 0.06;
  args_n2v.params[Neuron::g_ecs] = 0.55;
  Neuron n2v(args_n2v);
  n2v.set(Neuron::v, -67); n2v.set(Neuron::va, -67);
  n2v.set(Neuron::p, 1/(1+exp((-51-(-67.0))/10.3)));
  n2v.set(Neuron::q, 1/(1+exp((-45-(-67.0))/-3)));
  n2v.set(Neuron::h, 1/(1+exp((-55.2-(-67.0))/-7.1)));
  n2v.set(Neuron::n, 1.0/(1.0+exp((-30.0-(-67.0))/17.4)));

  Neuron::ConstructorArgs args_n3t;
  args_n3t.params[Neuron::n_type] = 3;
  args_n3t.params[Neuron::tau_p] = 4.0;
  args_n3t.params[Neuron::tau_q] = 400.0;
  args_n3t.params[Neuron::g_eca] = 8.0;
  args_n3t.params[Neuron::g_ecs] = 8.0;
  Neuron n3t(args_n3t);
  n3t.set(Neuron::v, -67.0); n3t.set(Neuron::va, -67.0);
  n3t.set(Neuron::p, 1/(1+exp((-61.6-(-67))/5.6)));
  n3t.set(Neuron::q, 1/(1+exp((-73.2-(-67))/-5.1)));
  n3t.set(Neuron::h, 1/(1+exp((-55.2-(-67.0))/-7.1)));
  n3t.set(Neuron::n, 1.0/(1.0+exp((-30.0-(-67.0))/17.4)));

  Neuron::ConstructorArgs args_so;
  args_so.params[Neuron::n_type] = 0;
  args_so.params[Neuron::tau_p] = 1.0;
  args_so.params[Neuron::tau_q] = 1.0;
  args_so.params[Neuron::g_eca] = 8.0;
  args_so.params[Neuron::g_ecs] = 8.0;
  Neuron so(args_so);
  so.set(Neuron::v, -67); so.set(Neuron::va, -67);
  so.set(Neuron::p, 0.0); so.set(Neuron::q, 0.0);
  so.set(Neuron::h, 1/(1+exp((-55.2-(-67.0))/-7.1)));
  so.set(Neuron::n, 1.0/(1.0+exp((-30.0-(-67.0))/17.4)));

  CGCNeuron::ConstructorArgs args_cgc;
  args_cgc.params[CGCNeuron::cm]=1.0;
  args_cgc.params[CGCNeuron::vna]=55.0; args_cgc.params[CGCNeuron::vk]=-90.0;
  args_cgc.params[CGCNeuron::vca]=80.0;
  args_cgc.params[CGCNeuron::Gnat]=1.68; args_cgc.params[CGCNeuron::Gnap]=0.44;
  args_cgc.params[CGCNeuron::Ga]=18.82; args_cgc.params[CGCNeuron::Gd]=1.20;
  args_cgc.params[CGCNeuron::Glva]=0.01; args_cgc.params[CGCNeuron::Ghva]=1.03;
  args_cgc.params[CGCNeuron::vh_h]=-56.43; args_cgc.params[CGCNeuron::vs_h]=-8.41;
  args_cgc.params[CGCNeuron::tau0_h]=778.82; args_cgc.params[CGCNeuron::delta_h]=0.03;
  args_cgc.params[CGCNeuron::vh_r]=-47.03; args_cgc.params[CGCNeuron::vs_r]=20.55;
  args_cgc.params[CGCNeuron::tau0_r]=4.01; args_cgc.params[CGCNeuron::delta_r]=1.00;
  args_cgc.params[CGCNeuron::vh_a]=-36.37; args_cgc.params[CGCNeuron::vs_a]=8.72;
  args_cgc.params[CGCNeuron::tau0_a]=13.28; args_cgc.params[CGCNeuron::delta_a]=0.39;
  args_cgc.params[CGCNeuron::vh_b]=-83.00; args_cgc.params[CGCNeuron::vs_b]=-6.20;
  args_cgc.params[CGCNeuron::tau0_b]=266.75; args_cgc.params[CGCNeuron::delta_b]=0.83;
  args_cgc.params[CGCNeuron::vh_n]=-59.43; args_cgc.params[CGCNeuron::vs_n]=34.79;
  args_cgc.params[CGCNeuron::tau0_n]=14.52; args_cgc.params[CGCNeuron::delta_n]=0.18;
  args_cgc.params[CGCNeuron::vh_e]=-14.25; args_cgc.params[CGCNeuron::vs_e]=6.96;
  args_cgc.params[CGCNeuron::tau0_e]=3.81; args_cgc.params[CGCNeuron::delta_e]=0.84;
  args_cgc.params[CGCNeuron::vh_f]=-21.44; args_cgc.params[CGCNeuron::vs_f]=-5.78;
  args_cgc.params[CGCNeuron::tau0_f]=34.68; args_cgc.params[CGCNeuron::delta_f]=0.97;
  args_cgc.params[CGCNeuron::Vh_m]=-35.20; args_cgc.params[CGCNeuron::Vs_m]=9.66;
  args_cgc.params[CGCNeuron::Vh_c]=-41.35; args_cgc.params[CGCNeuron::Vs_c]=5.05;
  args_cgc.params[CGCNeuron::Vh_d]=-64.13; args_cgc.params[CGCNeuron::Vs_d]=-4.03;
  CGCNeuron cgc(args_cgc);
  cgc.set(CGCNeuron::v, -60.0);
  cgc.set(CGCNeuron::h, 1.0/(1.0+exp((-56.43-(-60.0))/-8.41)));
  cgc.set(CGCNeuron::r, 1.0/(1.0+exp((-47.03-(-60.0))/20.55)));
  cgc.set(CGCNeuron::a, 1.0/(1.0+exp((-36.37-(-60.0))/8.72)));
  cgc.set(CGCNeuron::b, 1.0/(1.0+exp((-83.00-(-60.0))/-6.20)));
  cgc.set(CGCNeuron::n, 1.0/(1.0+exp((-59.43-(-60.0))/34.79)));
  cgc.set(CGCNeuron::e, 1.0/(1.0+exp((-14.25-(-60.0))/6.96)));
  cgc.set(CGCNeuron::f, 1.0/(1.0+exp((-21.44-(-60.0))/-5.78)));

  // =================== SINAPSIS CPG ===================
  // Sinapsis N1M -> N2v (EXCITATORIA LENTA)
  Synapse::ConstructorArgs syn_n1m_to_n2v;
  syn_n1m_to_n2v.params[Synapse::esyn] = 0.0;
  syn_n1m_to_n2v.params[Synapse::gsyn] = 0.077;
  syn_n1m_to_n2v.params[Synapse::tau_syn] = 200.0;
  syn_n1m_to_n2v.params[Synapse::v_pre] = -67.0;
  syn_n1m_to_n2v.params[Synapse::v_r] = -40.0;
  syn_n1m_to_n2v.params[Synapse::dec_slope] = 2.5;
  Synapse s_n1m_n2v(n1m, Neuron::v, n2v, Neuron::v, syn_n1m_to_n2v, 1);

  // Sinapsis N2v -> N1M (INHIBITORIA FUERTE)
  Synapse::ConstructorArgs syn_n2v_to_n1m;
  syn_n2v_to_n1m.params[Synapse::esyn] = -90.0;
  syn_n2v_to_n1m.params[Synapse::gsyn] = 50.0;
  syn_n2v_to_n1m.params[Synapse::tau_syn] = 50.0;
  syn_n2v_to_n1m.params[Synapse::v_pre] = -67.0;
  syn_n2v_to_n1m.params[Synapse::v_r] = -40.0;
  syn_n2v_to_n1m.params[Synapse::dec_slope] = 2.5;
  Synapse s_n2v_n1m(n2v, Neuron::v, n1m, Neuron::v, syn_n2v_to_n1m, 1);

  // Sinapsis N1M -> N3t (INHIBITORIA)
  Synapse::ConstructorArgs syn_n1m_to_n3t;
  syn_n1m_to_n3t.params[Synapse::esyn] = -90.0;
  syn_n1m_to_n3t.params[Synapse::gsyn] = 0.5;
  syn_n1m_to_n3t.params[Synapse::tau_syn] = 50.0;
  syn_n1m_to_n3t.params[Synapse::v_pre] = -67.0;
  syn_n1m_to_n3t.params[Synapse::v_r] = -40.0;
  syn_n1m_to_n3t.params[Synapse::dec_slope] = 2.5;
  Synapse s_n1m_n3t(n1m, Neuron::v, n3t, Neuron::v, syn_n1m_to_n3t, 1);

  // Sinapsis N3t -> N1M (INHIBITORIA)
  Synapse::ConstructorArgs syn_n3t_to_n1m;
  syn_n3t_to_n1m.params[Synapse::esyn] = -90.0;
  syn_n3t_to_n1m.params[Synapse::gsyn] = 8.0;
  syn_n3t_to_n1m.params[Synapse::tau_syn] = 50.0;
  syn_n3t_to_n1m.params[Synapse::v_pre] = -67.0;
  syn_n3t_to_n1m.params[Synapse::v_r] = -40.0;
  syn_n3t_to_n1m.params[Synapse::dec_slope] = 2.5;
  Synapse s_n3t_n1m(n3t, Neuron::v, n1m, Neuron::v, syn_n3t_to_n1m, 1);

  // Sinapsis N2v -> N3t (INHIBITORIA)
  Synapse::ConstructorArgs syn_n2v_to_n3t;
  syn_n2v_to_n3t.params[Synapse::esyn] = -90.0;
  syn_n2v_to_n3t.params[Synapse::gsyn] = 2.0;
  syn_n2v_to_n3t.params[Synapse::tau_syn] = 50.0;
  syn_n2v_to_n3t.params[Synapse::v_pre] = -67.0;
  syn_n2v_to_n3t.params[Synapse::v_r] = -40.0;
  syn_n2v_to_n3t.params[Synapse::dec_slope] = 2.5;
  Synapse s_n2v_n3t(n2v, Neuron::v, n3t, Neuron::v, syn_n2v_to_n3t, 1);

  // Sinapsis N2v -> SO (INHIBITORIA)
  Synapse::ConstructorArgs syn_n2v_to_so;
  syn_n2v_to_so.params[Synapse::esyn] = -90.0;
  syn_n2v_to_so.params[Synapse::gsyn] = 8.0;
  syn_n2v_to_so.params[Synapse::tau_syn] = 50.0;
  syn_n2v_to_so.params[Synapse::v_pre] = -67.0;
  syn_n2v_to_so.params[Synapse::v_r] = -40.0;
  syn_n2v_to_so.params[Synapse::dec_slope] = 2.5;
  Synapse s_n2v_so(n2v, Neuron::v, so, Neuron::v, syn_n2v_to_so, 1);

  // Sinapsis SO -> N1M (EXCITATORIA LENTA)
  Synapse::ConstructorArgs syn_so_to_n1m;
  syn_so_to_n1m.params[Synapse::esyn] = 0.0;
  syn_so_to_n1m.params[Synapse::gsyn] = 4.0;
  syn_so_to_n1m.params[Synapse::tau_syn] = 200.0;
  syn_so_to_n1m.params[Synapse::v_pre] = -67.0;
  syn_so_to_n1m.params[Synapse::v_r] = -40.0;
  syn_so_to_n1m.params[Synapse::dec_slope] = 2.5;
  Synapse s_so_n1m(so, Neuron::v, n1m, Neuron::v, syn_so_to_n1m, 1);

  // Sinapsis SO -> N2v (EXCITATORIA LENTA)
  Synapse::ConstructorArgs syn_so_to_n2v;
  syn_so_to_n2v.params[Synapse::esyn] = 0.0;
  syn_so_to_n2v.params[Synapse::gsyn] = 1.0;
  syn_so_to_n2v.params[Synapse::tau_syn] = 200.0;
  syn_so_to_n2v.params[Synapse::v_pre] = -67.0;
  syn_so_to_n2v.params[Synapse::v_r] = -40.0;
  syn_so_to_n2v.params[Synapse::dec_slope] = 2.5;
  Synapse s_so_n2v(so, Neuron::v, n2v, Neuron::v, syn_so_to_n2v, 1);

  // =================== SINAPSIS CGC ===================
  // Sinapsis CGC -> N1M (EXCITATORIA MONOSINÁPTICA)
  CGCSynapse::ConstructorArgs syn_cgc_to_n1m;
  syn_cgc_to_n1m.params[Synapse::esyn] = 0.0;
  syn_cgc_to_n1m.params[Synapse::gsyn] = 2.0;
  syn_cgc_to_n1m.params[Synapse::tau_syn] = 175.0;
  syn_cgc_to_n1m.params[Synapse::v_pre] = -67.0;
  syn_cgc_to_n1m.params[Synapse::v_r] = -40.0;
  syn_cgc_to_n1m.params[Synapse::dec_slope] = 2.5;
  CGCSynapse s_cgc_n1m(cgc, CGCNeuron::v, n1m, Neuron::v, syn_cgc_to_n1m, 1);

  // Sinapsis CGC -> N2v (EXCITATORIA POLISINÁPTICA)
  CGCSynapse::ConstructorArgs syn_cgc_to_n2v;
  syn_cgc_to_n2v.params[Synapse::esyn] = 0.0;
  syn_cgc_to_n2v.params[Synapse::gsyn] = 0.05;
  syn_cgc_to_n2v.params[Synapse::tau_syn] = 350.0;
  syn_cgc_to_n2v.params[Synapse::v_pre] = -67.0;
  syn_cgc_to_n2v.params[Synapse::v_r] = -40.0;
  syn_cgc_to_n2v.params[Synapse::dec_slope] = 2.5;
  CGCSynapse s_cgc_n2v(cgc, CGCNeuron::v, n2v, Neuron::v, syn_cgc_to_n2v, 1);

  // Sinapsis CGC -> SO (EXCITATORIA MONOSINÁPTICA)
  CGCSynapse::ConstructorArgs syn_cgc_to_so;
  syn_cgc_to_so.params[Synapse::esyn] = 0.0;
  syn_cgc_to_so.params[Synapse::gsyn] = 0.5;
  syn_cgc_to_so.params[Synapse::tau_syn] = 142.0;
  syn_cgc_to_so.params[Synapse::v_pre] = -67.0;
  syn_cgc_to_so.params[Synapse::v_r] = -40.0;
  syn_cgc_to_so.params[Synapse::dec_slope] = 2.5;
  CGCSynapse s_cgc_so(cgc, CGCNeuron::v, so, Neuron::v, syn_cgc_to_so, 1);

  // Sinapsis CGC -> N3t (EXCITATORIA MONOSINÁPTICA)
  CGCSynapse::ConstructorArgs syn_cgc_to_n3t;
  syn_cgc_to_n3t.params[Synapse::esyn] = 0.0;
  syn_cgc_to_n3t.params[Synapse::gsyn] = 0.05;
  syn_cgc_to_n3t.params[Synapse::tau_syn] = 350.0;
  syn_cgc_to_n3t.params[Synapse::v_pre] = -67.0;
  syn_cgc_to_n3t.params[Synapse::v_r] = -40.0;
  syn_cgc_to_n3t.params[Synapse::dec_slope] = 2.5;
  CGCSynapse s_cgc_n3t(cgc, CGCNeuron::v, n3t, Neuron::v, syn_cgc_to_n3t, 1);
  // =================== SIMULACIÓN ===================
  const double step = 0.01;
  const double t_transient = 5000.0;

  double half_ramp = n_ramp_steps * ramp_interval;
  double total_ramp = 2.0 * half_ramp;  // 1 triángulo completo
  const double simulation_time = t_transient + total_ramp + 5000.0;

  const int print_every = 10;
  int pcnt = 0;

  for (double time = 0; time < simulation_time; time += step) {
    s_n1m_n2v.step(step); s_n2v_n1m.step(step);
    s_n1m_n3t.step(step); s_n3t_n1m.step(step);
    s_n2v_n3t.step(step); s_n2v_so.step(step);
    s_so_n1m.step(step);  s_so_n2v.step(step);
    s_cgc_n1m.step(step); s_cgc_n2v.step(step);
    s_cgc_so.step(step);  s_cgc_n3t.step(step);

    // Calcular valor de la rampa (ya con signo correcto)
    double c_val = compute_ramp(time, t_transient, c_min, c_max,
                                n_ramp_steps, ramp_interval);

    // Inyectar corrientes directamente (sin negar, igual que cgc_cpg.cpp)
    so.add_synaptic_input( (ramp_target == 0) ? c_val : I_so);
    n1m.add_synaptic_input((ramp_target == 1) ? c_val : I_n1m);
    n2v.add_synaptic_input(I_n2v);
    n3t.add_synaptic_input((ramp_target == 3) ? c_val : I_n3t);
    cgc.add_synaptic_input((ramp_target == 4) ? c_val : I_cgc);

    n1m.step(step); n2v.step(step); n3t.step(step);
    so.step(step);  cgc.step(step);

    pcnt++;
    if (pcnt >= print_every) {
      pcnt = 0;
      std::cout << time << " "
                << n1m.get(Neuron::v) << " "
                << n2v.get(Neuron::v) << " "
                << n3t.get(Neuron::v) << " "
                << so.get(Neuron::v) << " "
                << cgc.get(CGCNeuron::v) << " "
                << c_val << std::endl;
    }
  }
  return 0;
}