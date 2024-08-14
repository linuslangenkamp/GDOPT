set(LIB_SRCS
        gdop.cpp
        integrator.cpp
        mesh.cpp
        problem.cpp
        solver.cpp
)

add_library(ipopt_do SHARED ${LIB_SRCS})
target_include_directories(ipopt_do PUBLIC ${CMAKE_BINARY_DIR}/third-party-install/usr/local/include/coin-or)
target_link_directories(ipopt_do PRIVATE ${CMAKE_BINARY_DIR}/third-party-install/usr/local/lib)
target_link_libraries(ipopt_do PRIVATE ipopt)
