cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
set(CMAKE_CXX_FLAGS "-O2")
project (acolyte)

find_package(Boost COMPONENTS system program_options filesystem timer REQUIRED)

# necessary for finding boost libraries on iridis 
include_directories(/local/software/boost/1.61.0/include)
find_library(boost_program_options boost_program_options /local/software/boost/1.61.0/lib)
find_library(boost_system boost_system /local/software/boost/1.61.0/lib)
find_library(boost_filesystem boost_filesystem /local/software/boost/1.61.0/lib)
find_library(boost_timer boost_filesystem /local/software/boost/1.61.0/lib)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

# find all source files
file(GLOB_RECURSE include include/*.hpp)
file(GLOB_RECURSE src src/*.cpp)
include_directories(include)
# include_directories(../delphes/install/include)
# find_library(Delphes_LIBRARIES Delphes ../delphes/install/lib)
include_directories(../delphes)
include_directories(../delphes/external)
find_library(Delphes_LIBRARIES Delphes ../delphes)

add_executable(analysis ${src})
target_link_libraries(analysis ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Delphes_LIBRARIES} ${boost_program_options} ${boost_system} ${boost_filesystem})
