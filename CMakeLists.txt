cmake_minimum_required(VERSION 3.22)
project(ipopt_do)

set(CMAKE_CXX_STANDARD 17)

include(ExternalProject)

ExternalProject_Add(mumps-ext
        GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Mumps.git
        GIT_TAG releases/3.0.7
        PATCH_COMMAND ${CMAKE_BINARY_DIR}/mumps-ext-prefix/src/mumps-ext/get.Mumps
        CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/third-party-install/usr/local
        BUILD_COMMAND make -j6
        BUILD_IN_SOURCE true)

ExternalProject_Add(ipopt-ext
        GIT_REPOSITORY https://github.com/coin-or/Ipopt.git
        GIT_TAG releases/3.14.16
        CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/third-party-install/usr/local
        BUILD_COMMAND make -j6
        BUILD_IN_SOURCE true
        DEPENDS mumps-ext)

file(GLOB SRC_FILES "src/*.cpp")

set(EXAMPLE_FILES "examples/hypersensitive.cpp"
        examples/hypersensitive.h)

include_directories(src examples)

add_executable(ipopt_do ${SRC_FILES} ${EXAMPLE_FILES})
target_include_directories(ipopt_do PRIVATE ${CMAKE_BINARY_DIR}/third-party-install/usr/local/include/coin-or)
target_link_directories(ipopt_do PRIVATE ${CMAKE_BINARY_DIR}/third-party-install/usr/local/lib)
target_link_libraries(ipopt_do ipopt)
get_property(dirs DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
message("Include directories: ${dirs}")
