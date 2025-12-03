/*************************************************************

Copyright (c) 2025, Alicia Garrido Peña <alicia.garrido@uam.es>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
    * Neither the name of the author nor the names of his contributors
      may be used to endorse or promote products derived from this
      software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*************************************************************/

#include<vector>
#include<string>
#include <DifferentialNeuronWrapper.h>
#include <ChemicalSynapsis.h>
#include <HindmarshRoseModel.h>
#include <SystemWrapper.h>
#include <RungeKutta4.h>
#include <iostream>
#include <cstring>
#include <map>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <filesystem>


using namespace std;

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<HindmarshRoseModel<double>>, Integrator> HR;


int main(int argc, char **argv) {
  if (argc < 3) {
      std::cerr << "Usage: " << argv[0]
                << " <output_file> <simulation_time> [step]\n";
      return 1;
  }

  // Required arguments
  std::string output_file = argv[1];
  double simulation_time  = std::atof(argv[2]);

  // Optional argument (step)
  double step = 0.01; // default value
  if (argc >= 4) {
      step = std::atof(argv[3]);
  }

  // Echo input back
  std::cout << "Output file: " << output_file << "\n";
  std::cout << "Simulation time: " << simulation_time << "\n";
  std::cout << "Step: " << step << "\n";

  // Initialize neuron model parameters
  HR::ConstructorArgs args;

  args.params[HR::e] = 0; // external current
  args.params[HR::mu] = 0.006;
  args.params[HR::S] = 4;
  args.params[HR::a] = 1;
  args.params[HR::b] = 3;
  args.params[HR::c] = 1;
  args.params[HR::d] = 5;
  args.params[HR::xr] = -1.6;
  args.params[HR::vh] = 1; // parameter for chaotic hyperpolarization

  HR h1(args);
  h1.set(HR::x, -0.712841);
  h1.set(HR::y, -1.93688);
  h1.set(HR::z, 3.16568);


  // Open output file
  std::ofstream out(output_file);
  if (!out.is_open()) {
      std::cerr << "Error: could not open " << output_file << " for writing.\n";
      return 1;
  }

  // Write header
  out << "Time V\n";

  // Simulation loop
  for (double time = 0; time < simulation_time; time += step) {
      h1.add_synaptic_input(2.5);
      h1.step(step);

      out << time << " "
          << h1.get(HR::x) << " "
          << "\n";
  }
  

  out.close();
  std::cout << "Simulation finished. Results written to " << output_file << "\n";


  return 0;
}
