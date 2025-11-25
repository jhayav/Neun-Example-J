cmake_minimum_required(VERSION 3.5)

project(HR LANGUAGES CXX)

# Configuración del estándar de C++
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Flags similares al Makefile
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -I/usr/local/Neun/0.4.0")

# Ejecutable
add_executable(HR
    HR.cpp
)

# Enlazado equivalente a LDLIBS = -lyaml-cpp
target_link_libraries(HR
    yaml-cpp
)
