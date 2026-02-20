#include <DifferentialNeuronWrapper.h>
#include <VavoulisModel.h>
#include <RungeKutta4.h>
#include <SystemWrapper.h>
#include <iostream>

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<VavoulisModel<double>>, Integrator>
    Neuron;

int main(int argc, char **argv) {

  Neuron::ConstructorArgs args;

  args.params[Neuron::n_type] = 2;     // n2v
  // Para N2v, tau_p y tau_q se calculan dinamicamente en el modelo
  // segun las ecuaciones de la Tabla 1:
  // tau_p = 28.3 + 44.1 * exp(-((-11.8-V_A)/26.6)^2)
  // tau_q = 187.6 + 637.7 * exp(-((-9.5-V_A)/23.3)^2)
  // Estos valores de parametro no se usan para N2v (se ignoran)
  args.params[Neuron::tau_p] = 1.0;    // No usado para N2v
  args.params[Neuron::tau_q] = 1.0;    // No usado para N2v
  args.params[Neuron::g_eca] = 0.06;   // Acoplamiento axon->soma (muy debil)
  args.params[Neuron::g_ecs] = 0.55;   // Acoplamiento soma->axon

  // Initialize a new neuron model
  Neuron n(args);

  n.set(Neuron::v, -67);
  n.set(Neuron::va, -67);
  n.set(Neuron::p, 1 / (1 + exp((-51 - (-67.0))/10.3)));
  n.set(Neuron::q, 1 / (1 + exp((-45 - (-67.0))/-3)));
  n.set(Neuron::h, 1 / (1 + exp((-55.2 - (-67.0))/-7.1)));
  n.set(Neuron::n, 1.0 / (1.0 + exp((-30.0 - (-67.0)) / 17.4)));

  // Parametros de simulacion
  const double step = 0.01;            // Paso de integracion (ms)
  const double simulation_time = 3000; // Tiempo total (ms) - 3 segundos
  
  const double t_pulse_start = 300;    // Inicio del pulso (ms)
  const double t_pulse_end = 2700;     // Fin del pulso (ms)
  const double I_inj = -5.0;           // Corriente despolarizante

  // Perform the simulation
  for (double time = 0; time < simulation_time; time += step) {
    if (time >= t_pulse_start && time <= t_pulse_end) {
      n.add_synaptic_input(I_inj);
    }
    
    n.step(step);
    // Salida: tiempo, V_soma, V_axon, p, q, h, n
    std::cout << time << " " << n.get(Neuron::v) << " " << n.get(Neuron::va) 
              << " " << n.get(Neuron::p) << " " << n.get(Neuron::q)
              << " " << n.get(Neuron::h) << " " << n.get(Neuron::n) << std::endl;
  }

  return 0;
}
