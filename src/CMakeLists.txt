set(LIB_SRCS
        gdop.cpp
        integrator.cpp
        mesh.cpp
        problem.cpp
        solver.cpp
)

add_library(gdopt SHARED ${LIB_SRCS})
target_include_directories(gdopt PRIVATE ${CMAKE_BINARY_DIR}/third-party-install/usr/local/include/coin-or)
target_link_directories(gdopt PRIVATE ${CMAKE_BINARY_DIR}/third-party-install/usr/local/lib)
target_link_libraries(gdopt PRIVATE ipopt coinmumps)
target_link_libraries(gdopt PRIVATE dl)
target_link_libraries(gdopt PRIVATE lapack)
target_compile_options(gdopt PRIVATE -ffast-math)
add_dependencies(gdopt ipopt-ext)


set(GDOPT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})

if (PY_BUILD_CMAKE_MODULE_NAME)
        set(GDOPT_INSTALL_PREFIX ${PY_BUILD_CMAKE_MODULE_NAME})
endif ()

install(TARGETS gdopt 
        LIBRARY DESTINATION ${GDOPT_INSTALL_PREFIX}/lib
        COMPONENT libgdopt)

install(FILES constants.h integrator.h mesh.h problem.h gdop.h solver.h constraint.h expression.h
        DESTINATION ${GDOPT_INSTALL_PREFIX}/include/libgdopt/
        COMPONENT libgdopt)