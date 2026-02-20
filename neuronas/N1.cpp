#include <DifferentialNeuronWrapper.h>
#include <VavoulisModel.h>
#include <RungeKutta4.h>
#include <SystemWrapper.h>
#include <CurrentPulse.h>
#include <iostream>

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<VavoulisModel<double>>, Integrator>
    Neuron;

int main(int argc, char **argv) {
  // Struct to initialize neuron model parameters
  Neuron::ConstructorArgs args;

  // Set the parameter values (Tabla 1, Vavoulis 2007)
  args.params[Neuron::n_type] = 1; // n1m
  args.params[Neuron::tau_p] = 250.0;  // Constante de tiempo lenta (ms)
  args.params[Neuron::tau_q] = 1.0;    // No se usa para N1M
  args.params[Neuron::g_eca] = 8.0;    // Acoplamiento axon->soma
  args.params[Neuron::g_ecs] = 8.0;    // Acoplamiento soma->axon

  // Initialize a new neuron model
  Neuron n(args);

  // Condiciones iniciales
  n.set(Neuron::v, -67);
  n.set(Neuron::va, -67);
  n.set(Neuron::p, 1 / (1 + exp((-38.8 - (-67.0)) / 10.0)));
  n.set(Neuron::q, 0.0);
  n.set(Neuron::h, 1 / (1 + exp((-55.2 - (-67.0))/-7.1)));
  n.set(Neuron::n, 1.0 / (1.0 + exp((-30.0 - (-67.0)) / 17.4)));

  // Parametros de simulacion
  const double step = 0.01;           // Paso de integracion (ms)
  const double simulation_time = 2000; // Tiempo total (ms) - 2 segundos
  
  // Parametros de estimulacion (Fig. 3A del paper)
  // Pulso de corriente despolarizante para activar el plateau
  // En el modelo, SYNAPTIC_INPUT tiene signo negativo en la ecuacion,
  // por lo que una corriente NEGATIVA produce despolarizacion
  const double t_pulse_start = 200;   // Inicio del pulso (ms)
  const double t_pulse_end = 1800;    // Fin del pulso (ms)
  const double I_inj = -10.0;         // Corriente despolarizante moderada

  // Perform the simulation
  for (double time = 0; time < simulation_time; time += step) {
    // Inyeccion de corriente durante el pulso
    if (time >= t_pulse_start && time <= t_pulse_end) {
      n.add_synaptic_input(I_inj);
    }
    
    n.step(step);
    // Salida: tiempo, V_soma, V_axon, p, h, n
    std::cout << time << " " << n.get(Neuron::v) << " " << n.get(Neuron::va) 
              << " " << n.get(Neuron::p) << " " << n.get(Neuron::h) 
              << " " << n.get(Neuron::n) << std::endl;
  }

  return 0;
}
