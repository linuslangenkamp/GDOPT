cmake_minimum_required(VERSION 3.22)
project(ipopt_do)

set(CMAKE_CXX_STANDARD 17)

add_executable(ipopt_do src/main.cpp
        src/test_nlp.cpp
        src/test_nlp.h
        src/integrator.cpp
        src/integrator.h
        src/problem.h
        src/expression.h
        src/constraint.h
        src/mesh.h
        src/gdop.h
        src/gdop.cpp
        src/mesh.cpp
        src/constants.h
        src/problem.cpp)
target_include_directories(ipopt_do PRIVATE /usr/local/include/coin-or)
target_link_directories(ipopt_do PRIVATE /usr/local/lib)
target_link_libraries(ipopt_do ipopt coinmumps)