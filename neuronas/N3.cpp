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

  args.params[Neuron::n_type]  = 3;     // n3t
  
  // Constantes de tiempo para N3t (Tabla 1)
  // tau_p = 4 ms (muy rapida, activacion del canal T)
  // tau_q = 400 ms (lenta, inactivacion del canal T)
  args.params[Neuron::tau_p] = 4.0;    // ms
  args.params[Neuron::tau_q] = 400.0;  // ms
  args.params[Neuron::g_eca] = 8.0;    // Acoplamiento axon->soma
  args.params[Neuron::g_ecs] = 8.0;    // Acoplamiento soma->axon

  Neuron n(args);

  n.set(Neuron::v, -67.0);   // Soma en potencial de fuga
  n.set(Neuron::va, -67.0);  // Axón en potencial de fuga
  // p_inf = 1 / (1 + exp((-61.6 - (-67))/5.6)) = 0.276
  n.set(Neuron::p, 1 / (1 + exp((-61.6 - (-67))/5.6)));   
  // q_inf = 1 / (1 + exp((-73.2 - (-67))/-5.1)) = 0.229
  n.set(Neuron::q, 1 / (1 + exp((-73.2 - (-67))/-5.1)));   
  n.set(Neuron::h, 1 / (1 + exp((-55.2 - (-67.0))/-7.1)));
  n.set(Neuron::n, 1.0 / (1.0 + exp((-30.0 - (-67.0)) / 17.4)));

  // Parametros de simulacion
  const double step = 0.01;            // Paso de integracion (ms)
  const double simulation_time = 4000; // Tiempo total (ms) - 4 segundos
  
  // Parametros de estimulacion (Fig. 3C del paper)
  // N3t muestra post-inhibitory rebound con pulsos hiperpolarizantes breves
  // Segun Fig. 3C(iii): iinj = -2, -4, -8 mV
  const double t_pulse_start = 1000;
  const double t_pulse_end = 1800;
  const double I_inj = 8.0;            // Corriente hiperpolarizante (positiva = hiperpolariza en el modelo)

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
