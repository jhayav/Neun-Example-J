set (INCLUDE_DIR ../include)
include_directories(${INCLUDE_DIR} ../concepts ../models ../integrators ../wrappers ../archetypes)

add_executable(N1M N1.cpp)
target_link_libraries(N1M)

add_executable(N2v N2.cpp)
target_link_libraries(N2v)

add_executable(N3t N3.cpp)
target_link_libraries(N3t)

add_executable(SO SO.cpp)
target_link_libraries(SO)

add_executable(CGC CGC.cpp)
target_link_libraries(CGC)
