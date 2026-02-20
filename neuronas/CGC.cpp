/*************************************************************
 * CGC.cpp - Neurona CGC (Cerebral Giant Cell) aislada
 * 
 * Modelo de un compartimento
 * con 6 corrientes iónicas:
 *   - iNaT: Sodio transitorio rápido (Gnat, h)
 *   - iNaP: Sodio persistente (Gnap, r)
 *   - iA:   Potasio tipo A (Ga, a, b)
 *   - iD:   Potasio rectificador retardado (Gd, n)
 *   - iLVA: Calcio de bajo umbral (Glva, c_inf, d_inf)
 *   - iHVA: Calcio de alto umbral (Ghva, e, f)
 * 
 * Variables dinámicas: v, h, r, a, b, n, e, f
 * 
 *************************************************************/

#include <DifferentialNeuronWrapper.h>
#include <VavoulisCGCModel.h>
#include <RungeKutta4.h>
#include <SystemWrapper.h>
#include <iostream>

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<VavoulisCGCModel<double>>, Integrator>
    Neuron;

int main(int argc, char **argv) {

  Neuron::ConstructorArgs args;

  // PARÁMETROS DEL MODELO (Tabla 1)

  // Capacitancia y potenciales de reversión (mV)
  args.params[Neuron::cm]  = 1.0;      // Capacitancia de membrana (µF/cm²)
  args.params[Neuron::vna] = 55.0;     // E_Na (mV)
  args.params[Neuron::vk]  = -90.0;    // E_K (mV)
  args.params[Neuron::vca] = 80.0;     // E_Ca (mV)

  // Conductancias máximas (mS/cm²)
  args.params[Neuron::Gnat] = 1.68;    // Sodio transitorio
  args.params[Neuron::Gnap] = 0.44;    // Sodio persistente
  args.params[Neuron::Ga]   = 18.82;   // Potasio tipo A
  args.params[Neuron::Gd]   = 1.20;    // Potasio rectificador retardado
  args.params[Neuron::Glva] = 0.01;    // Calcio LVA (bajo umbral)
  args.params[Neuron::Ghva] = 1.03;    // Calcio HVA (alto umbral)

  //  h: inactivación Na transitorio 
  args.params[Neuron::vh_h]    = -56.43;  // V_half (mV)
  args.params[Neuron::vs_h]    = -8.41;   // Pendiente (mV), negativa = inactivación
  args.params[Neuron::tau0_h]  = 778.82;  // Constante de tiempo basal (ms)
  args.params[Neuron::delta_h] = 0.03;    // Factor de asimetría de tau (adim.)

  //  r: activación Na persistente 
  args.params[Neuron::vh_r]    = -47.03;  // V_half (mV)
  args.params[Neuron::vs_r]    = 20.55;   // Pendiente (mV), positiva = activación
  args.params[Neuron::tau0_r]  = 4.01;    // Constante de tiempo basal (ms)
  args.params[Neuron::delta_r] = 1.00;    // Factor de asimetría de tau (adim.)

  //  a: activación K tipo A 
  args.params[Neuron::vh_a]    = -36.37;  // V_half (mV)
  args.params[Neuron::vs_a]    = 8.72;    // Pendiente (mV), positiva = activación
  args.params[Neuron::tau0_a]  = 13.28;   // Constante de tiempo basal (ms)
  args.params[Neuron::delta_a] = 0.39;    // Factor de asimetría de tau (adim.)

  //  b: inactivación K tipo A 
  args.params[Neuron::vh_b]    = -83.00;  // V_half (mV)
  args.params[Neuron::vs_b]    = -6.20;   // Pendiente (mV), negativa = inactivación
  args.params[Neuron::tau0_b]  = 266.75;  // Constante de tiempo basal (ms)
  args.params[Neuron::delta_b] = 0.83;    // Factor de asimetría de tau (adim.)

  //  n: activación K rectificador retardado 
  args.params[Neuron::vh_n]    = -59.43;  // V_half (mV)
  args.params[Neuron::vs_n]    = 34.79;   // Pendiente (mV), positiva = activación
  args.params[Neuron::tau0_n]  = 14.52;   // Constante de tiempo basal (ms)
  args.params[Neuron::delta_n] = 0.18;    // Factor de asimetría de tau (adim.)

  //  e: activación Ca alto umbral (HVA) 
  args.params[Neuron::vh_e]    = -14.25;  // V_half (mV)
  args.params[Neuron::vs_e]    = 6.96;    // Pendiente (mV), positiva = activación
  args.params[Neuron::tau0_e]  = 3.81;    // Constante de tiempo basal (ms)
  args.params[Neuron::delta_e] = 0.84;    // Factor de asimetría de tau (adim.)

  //  f: inactivación Ca alto umbral (HVA) 
  args.params[Neuron::vh_f]    = -21.44;  // V_half (mV)
  args.params[Neuron::vs_f]    = -5.78;   // Pendiente (mV), negativa = inactivación
  args.params[Neuron::tau0_f]  = 34.68;   // Constante de tiempo basal (ms)
  args.params[Neuron::delta_f] = 0.97;    // Factor de asimetría de tau (adim.)

  //  m: activación Na transitorio (instantánea, sin dinámica propia) 
  args.params[Neuron::Vh_m] = -35.20;    // V_half (mV)
  args.params[Neuron::Vs_m] = 9.66;      // Pendiente (mV)

  //  c, d: activación/inactivación Ca bajo umbral LVA (instantáneas) 
  args.params[Neuron::Vh_c] = -41.35;    // V_half activación (mV)
  args.params[Neuron::Vs_c] = 5.05;      // Pendiente activación (mV)
  args.params[Neuron::Vh_d] = -64.13;    // V_half inactivación (mV)
  args.params[Neuron::Vs_d] = -4.03;     // Pendiente inactivación (mV)

  // Inicializar neurona
  Neuron n(args);


  // Condiciones iniciales (reposo basal aprox -60 mV)
  // x_inf(V) = 1/(1+exp((vh-V)/vs)), evaluadas en V = -60 mV
  n.set(Neuron::v, -60.0);
  n.set(Neuron::h, 1.0 / (1.0 + exp((-56.43 - (-60.0)) / -8.41)));
  n.set(Neuron::r, 1.0 / (1.0 + exp((-47.03 - (-60.0)) / 20.55)));
  n.set(Neuron::a, 1.0 / (1.0 + exp((-36.37 - (-60.0)) / 8.72)));
  n.set(Neuron::b, 1.0 / (1.0 + exp((-83.00 - (-60.0)) / -6.20)));
  n.set(Neuron::n, 1.0 / (1.0 + exp((-59.43 - (-60.0)) / 34.79)));
  n.set(Neuron::e, 1.0 / (1.0 + exp((-14.25 - (-60.0)) / 6.96)));
  n.set(Neuron::f, 1.0 / (1.0 + exp((-21.44 - (-60.0)) / -5.78)));


  const double step = 0.01;              // Paso de integración (ms)
  const double simulation_time = 3000;   // Tiempo total (ms)

  // Estimulación (rango típico: 0-2 nA)
  const double t_pulse_start = 500.0;    // Inicio del pulso (ms)
  const double t_pulse_end   = 2500.0;   // Fin del pulso (ms)
  const double I_inj         = 0.2;     // Corriente inyectada

  // Simulación
  for (double time = 0; time < simulation_time; time += step) {
    if (time >= t_pulse_start && time <= t_pulse_end) {
      n.add_synaptic_input(I_inj);
    }

    n.step(step);

    // Salida: tiempo, V, h, r, a, b, n, e, f
    std::cout << time << " "
              << n.get(Neuron::v) << " "
              << n.get(Neuron::h) << " "
              << n.get(Neuron::r) << " "
              << n.get(Neuron::a) << " "
              << n.get(Neuron::b) << " "
              << n.get(Neuron::n) << " "
              << n.get(Neuron::e) << " "
              << n.get(Neuron::f)
              << std::endl;
  }

  return 0;
}