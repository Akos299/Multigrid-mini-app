# file( GLOB HEADER_LIST ../include/*.hpp)
# file(GLOB SOURCE_LIST multigrid_base.cpp mgrid_stack)
set(HEADER_LIST ../include/utilities.hpp
../include/transfert_operators.hpp
../include/multigrid_base.hpp
../include/mgrid_stack.hpp
../include/boundary_condition.hpp)

set(SOURCE_LIST multigrid_base.cpp mgrid_stack.cpp)
add_library(multigrid_library ${SOURCE_LIST} ${HEADER_LIST})
target_include_directories(multigrid_library PUBLIC ../include)
target_compile_features(multigrid_library PUBLIC cxx_std_17)
set_target_properties(multigrid_library
    PROPERTIES COMPILER_FLAGS "-g -m64 -arch x86_64 -msse -Wall -pedantic"
               LINKER_FLAGS "-g -m64 -arch x86_64 -msse")
