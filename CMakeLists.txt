cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project (apsis)

find_package(Boost COMPONENTS system program_options filesystem REQUIRED)
if(NOT Boost_FOUND)
    include_directories(/local/software/boost/1.61.0/include)
    find_library(lboost_program_options lboost_program_options /local/software/boost/1.61.0/lib)
    find_library(lboost_system lboost_system /local/software/boost/1.61.0/lib)
    find_library(lboost_filesystem lboost_filesystem /local/software/boost/1.61.0/lib)
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

# find all source files
file(GLOB_RECURSE include include/*.hpp)
file(GLOB_RECURSE src src/*.cpp)
include_directories(include)
include_directories(../delphes/install/include)
find_library(Delphes_LIBRARIES Delphes ../delphes/install/lib)

add_executable(analysis ${src})
target_link_libraries(analysis ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Delphes_LIBRARIES} ${lboost_program_options} ${lboost_system} ${lboost_filesystem})
