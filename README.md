# Neun-Example
Example of a modelling project using [Neun](github.com/GNB-UAM/Neun)
This is a repository with several examples of C++ code using models available in Neun library. 
## 0. Install Neun
Neun is a library, so in order to use it in your project you can: 

  a. (Recommended) Install it in your Linux system following the instructions at github.com/GNB-UAM/Neun
  
  b. Clone it and reference the corresponding folder file
  
## 1. Create your Makefile 
Since this project is in C++ we need to compile it. For that purpose, CMake helps us create a Makefile that will generate the necessary contents.

    mkdir build; cd build 
    cmake ..

Note that now you have a complete Makefile including all the examples. 

You can now run ``make`` to generate the executable files
   
    make 

After make you should have a new file HR. You can try running it doing

    ./HR 

Include arguments to generate your first model, for example:

    ./HR test.dat 10000 0.01

And you can visualize it by running:

    python ../plot.py test.dat

Note that you may need to install some common package for using python. You can create an environment by doing this:
  
    python -m venv neun-py-env
    source neun-py-env/bin/activate
    pip install -r ../requirements.txt

## Including your own cpp [Coming soon]
Use template.cpp to generate the model you prefer. You will need to choose:

  1. Model: changing this include in the wrapper
  2. Integrator:
  3. Synapse (In case you want to simulate a Network)

Once your includes are ready tune all the necessary parameters (check the .h file of your corresponding model)

In the template.cpp there is a for loop that performs the following steps:

  1. includes an external input in case it is neccessary neuron.input(1.5)
  2. updates the neural model by performing neuron.step()
  3. Prints the updated V value using the default output channel (cout)

Note that when running this file you will have all the values of your simulation and you can write that into a new file to plot it in your preferred tool.

plot.py is an example of code to plot your simulation from a file using Python. 

