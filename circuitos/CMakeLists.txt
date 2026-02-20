set (INCLUDE_DIR ../include)
include_directories(${INCLUDE_DIR} ../concepts ../models ../integrators ../wrappers ../archetypes)

add_executable(N1N2 n1-2.cpp)
target_link_libraries(N1N2)

add_executable(N1N2N3 n1-2-3.cpp)
target_link_libraries(N1N2N3)

add_executable(CPG cpg_completo.cpp)
target_link_libraries(CPG)