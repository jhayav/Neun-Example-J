set (INCLUDE_DIR ../include)
include_directories(${INCLUDE_DIR} ../concepts ../models ../integrators ../wrappers ../archetypes)

add_executable(N1N2 n1-2.cpp)
target_link_libraries(N1N2)

add_executable(N1N2N3 n1-2-3.cpp)
target_link_libraries(N1N2N3)

add_executable(CPG cpg_completo.cpp)
target_link_libraries(CPG)

add_executable(cgc cgc_cpg.cpp)
target_link_libraries(cgc)

add_executable(circuito_gen circuito_gen.cpp)
target_link_libraries(circuito_gen)

# add_executable(cgc_valores cgc_cpg_valores.cpp)
# target_link_libraries(cgc_valores)