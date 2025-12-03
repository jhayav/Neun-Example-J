cmake_minimum_required(VERSION 3.5)

project(HR LANGUAGES CXX)

# Configuration of the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Flags for the compiler
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -I/usr/local/Neun/0.4.0")

# Executable
add_executable(HR HR.cpp)
add_executable(basic basic.cpp)
add_executable(synapsis synapsis.cpp)
add_executable(chemicalSynapsis chemicalSynapsis.cpp)

# Linking equivalent to LDLIBS = -lyaml-cpp
# target_link_libraries(HR
#     yaml-cpp
# )
