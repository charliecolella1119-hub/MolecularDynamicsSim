cmake_minimum_required(VERSION 3.16)

project(MolecularDynamicsSim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SFML 3 COMPONENTS Graphics Window System REQUIRED)

add_executable(MolecularDynamicsSim
    src/main.cpp
)

target_link_libraries(MolecularDynamicsSim
    SFML::Graphics
    SFML::Window
    SFML::System
)