add_executable(poisson poisson.cpp main.cpp poisson.hpp)

target_link_libraries(poisson PUBLIC multigrid_library)
target_compile_features(poisson PUBLIC cxx_std_17)
target_include_directories(poisson PUBLIC ../include )
# set_target_properties(poisson
    # PROPERTIES COMPILER_FLAGS "-g -m64 -arch x86_64 -msse -Wall -pedantic -v"
    #            LINKER_FLAGS "-g -m64 -arch x86_64 -msse -v")