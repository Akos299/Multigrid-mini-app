# set the headers list
set(HEADER_LIST ../include/utilities.hpp
../include/transfert_operators.hpp
../include/multigrid_base.hpp
../include/mgrid_stack.hpp
../include/boundary_condition.hpp
../include/ndArray.h)


# set the sources
set(SOURCE_LIST multigrid_base.cpp mgrid_stack.cpp)
#
add_library(multigrid_library ${SOURCE_LIST})
# point target to headers directory
target_include_directories(multigrid_library PUBLIC ../include)
#set C++ version
target_compile_features(multigrid_library PUBLIC cxx_std_17)
# set compiler options
set_target_properties(multigrid_library
    PROPERTIES COMPILER_FLAGS "-g -m64 -arch x86_64 -msse -Wall -pedantic"
               LINKER_FLAGS "-g -m64 -arch x86_64 -msse")
