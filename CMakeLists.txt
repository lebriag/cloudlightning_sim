cmake_minimum_required(VERSION 2.8.1)
project(CloudLightning)

###############################################################################
## file globbing ##############################################################
###############################################################################

# These instructions search the directory tree when cmake is
# invoked and put all files that match the pattern in the variable
# 'sources'.
file(GLOB_RECURSE sources      src/*.cpp include/*.h)
# You can use set(sources src/main.cpp) etc if you don't want to
# use globbing to find files automatically.

find_package(MPI REQUIRED)
if (MPI_FOUND)
  SET(CMAKE_C_COMPILER mpicc)
  SET(CMAKE_CXX_COMPILER mpicxx)
endif()

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory(doc)

###############################################################################
## target definitions #########################################################
###############################################################################

# The data is just added to the executable, because in some IDEs (QtCreator)
# files are invisible when they are not explicitly part of the project.
add_executable(cl_sim ${sources} ${data})

# Just for example add some compiler flags.
target_compile_options(cl_sim PUBLIC -O3 -std=c++11 -lm -Wall)

# This allows to include files relative to the root of the src directory with a <> pair
target_include_directories(cl_sim PUBLIC include/)

target_link_libraries(cl_sim ${MPI_C_LIBRARIES})

# Added support for the include what you use tool
find_program(iwyu_path NAMES include-what-you-use iwyu)
if(iwyu_path)
  set_property(TARGET cl_sim PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path})
endif()

install (TARGETS cl_sim DESTINATION ${CMAKE_BINARY_DIR}/bin)
