set (INCLUDE_DIR ../include)
include_directories(${INCLUDE_DIR} ../concepts ../models ../integrators ../wrappers ../archetypes)

add_executable(basic basic.cpp)
add_executable(synapsis synapsis.cpp)
