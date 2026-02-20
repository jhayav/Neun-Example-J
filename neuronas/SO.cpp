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

  args.params[Neuron::n_type] = 0;
  
  // SO es pasiva, no tiene corrientes intrinsecas (ix = 0)
  // tau_p y tau_q no se usan para SO
  args.params[Neuron::tau_p] = 1.0;    // No usado
  args.params[Neuron::tau_q] = 1.0;    // No usado
  
  // Acoplamiento electrico (igual que N1M y N3t)
  // g_ec = 8 segun Tabla 1
  args.params[Neuron::g_eca] = 8.0;    // Acoplamiento axon->soma
  args.params[Neuron::g_ecs] = 8.0;    // Acoplamiento soma->axon

  Neuron n(args);

  // Condiciones iniciales
  n.set(Neuron::v, -67);
  n.set(Neuron::va, -67);
  n.set(Neuron::p, 0.0);    // No usado para SO
  n.set(Neuron::q, 0.0);    // No usado para SO
  n.set(Neuron::h, 1 / (1 + exp((-55.2 - (-67.0))/-7.1)));
  n.set(Neuron::n, 1.0 / (1.0 + exp((-30.0 - (-67.0)) / 17.4)));

  // Parametros de simulacion
  const double step = 0.01;            // Paso de integracion (ms)
  const double simulation_time = 3000; // Tiempo total (ms) - 3 segundos
  
  // Parametros de estimulacion (Fig. 3D del paper)
  // SO es pasiva: dispara con corriente despolarizante
  const double t_stim_start = 200;     // Inicio estimulacion continua
  const double t_stim_end = 2800;      // Fin estimulacion
  const double I_stim = -10.0;         // Corriente despolarizante
  
  // Pulsos inhibitorios breves (como en Fig. 3D(ii))
  const double t_inhib1_start = 800;
  const double t_inhib1_end = 900;
  const double t_inhib2_start = 1500;
  const double t_inhib2_end = 1600;
  const double t_inhib3_start = 2200;
  const double t_inhib3_end = 2300;
  const double I_inhib = 15.0;         // Corriente hiperpolarizante

  for (double time = 0; time < simulation_time; time += step) {
    double current = 0.0;
    
    // Estimulacion despolarizante continua
    if (time >= t_stim_start && time <= t_stim_end) {
      current = I_stim;
    }
    
    // Pulsos inhibitorios breves
    if ((time >= t_inhib1_start && time <= t_inhib1_end) ||
        (time >= t_inhib2_start && time <= t_inhib2_end) ||
        (time >= t_inhib3_start && time <= t_inhib3_end)) {
      current = I_inhib;
    }
    
    n.add_synaptic_input(current);
    n.step(step);
    
    // Salida: tiempo, V_soma, V_axon, p, q, h, n
    std::cout << time << " " << n.get(Neuron::v) << " " << n.get(Neuron::va) 
              << " " << n.get(Neuron::p) << " " << n.get(Neuron::q)
              << " " << n.get(Neuron::h) << " " << n.get(Neuron::n) << std::endl;
  }

  return 0;
}
