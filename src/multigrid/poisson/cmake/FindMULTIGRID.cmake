# Find multigrid include directories and libraries
#
# MULTIGRID_INCLUDE_DIRECTORIES - where to find netcdf.h
# MULTIGRID_LIBRARIES - list of libraries to link against when using NetCDF
# MULTIGRID_FOUND - Do not attempt to use NetCDF if "no", "0", or undefined.

include(LibFindMacros)

# Dependencies 
set(MULTIGRID_PREFIX "/usr/local" 
    CACHE PATH "Path to search for Multigrid header and library files" )
# ==================
# = Find Multigrid =
# ==================
# Include dir
find_path(MULTIGRID_INCLUDE_DIR NAMES multigrid/multigrid_base.hpp 
    PATHS ${MULTIGRID_PREFIX})

# Finally the library itself
find_library(MULTIGRID_LIBRARY NAMES multigrid PATHS ${MULTIGRID_PREFIX})

# Set the include dir variables and the libraries and let libfind_process do
# the rest. NOTE: Singular variables for this library, plural for libraries
# this this lib depends on.
set(MULTIGRID_PROCESS_INCLUDES 
    MULTIGRID_INCLUDE_DIR 
)
set(MULTIGRID_PROCESS_LIBS 
    MULTIGRID_LIBRARY 
    )
libfind_process(MULTIGRID)