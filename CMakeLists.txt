Project (PARADIS)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_VERBOSE_MAKEFILE OFF)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

message ("** PARADIS_BINARY_DIR =  ${PARADIS_BINARY_DIR}")
message ("** PARADIS_SOURCE_DIR =  ${PARADIS_SOURCE_DIR}")


#find_package(MPI REQUIRED)
find_package(MPI)
if(MPI_FOUND)  
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_EXE_LINKER_FLAGS}")  
endif()  

set(CMAKE_C_FLAGS_GDB "-ggdb" CACHE STRING "C flag for gdb" )
set(CMAKE_CXX_FLAGS_GDB "-ggdb" CACHE STRING "CXX flag for gdb" )
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." )
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "GDB")
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
message("** CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

message("** CMAKE_C_COMPILER = ${CMAKE_C_COMPILER} (${CMAKE_C_FLAGS})")
message("** CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_FLAGS}); CMAKE_EXE_LINKER_FLAGS = ${CMAKE_EXE_LINKER_FLAGS}")

Option(INSTALL_IN_PARADIS_SOURCE_DIR "install executable files in ${PARADIS_SOURCE_DIR}" ON)
if(INSTALL_IN_PARADIS_SOURCE_DIR)
    Set(CMAKE_INSTALL_PREFIX ${PARADIS_SOURCE_DIR} CACHE PATH "Install directory" FORCE)
endif(INSTALL_IN_PARADIS_SOURCE_DIR)
message("** CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")

set(PARADIS_EXE paradis CACHE STRING "ParaDiS executable file")
Set(PARADISGEN_EXE paradisgen CACHE STRING "ParaDiS generate executable file" )
message("** PARADIS_EXE = ${PARADIS_EXE}; PARADISGEN_EXE = ${PARADISGEN_EXE}")

Option(OPENMP_MODE "Build paradis using definition OPENMP_MODE" OFF)
if(OPENMP_MODE)
    find_package(OpenMP)  
    if(OPENMP_FOUND)  
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")  
        message("** Find OpenMP: OpenMP_C_FLAGS = ${OpenMP_C_FLAGS}; OpenMP_CXX_FLAGS =  ${OpenMP_CXX_FLAGS}; OpenMP_EXE_LINKER_FLAGS = ${OpenMP_EXE_LINKER_FLAGS}")
    else(OPENMP_FOUND)
        add_definitions(-DOPENMP_MODE=OFF)
    endif(OPENMP_FOUND)  
endif(OPENMP_MODE)

Find_Package(GEOS)
IF(GEOS_FOUND)
    message(STATUS "include and lib " ${GEOS_INCLUDE_DIR} " " ${GEOS_LIB_DIR})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GEOS_C_FLAGS}")  
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GEOS_CXX_FLAGS}")  
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GEOS_EXE_LINKER_FLAGS}")  
    message(STATUS "GEOS_C_FLAGS = " ${GEOS_C_FLAGS})  
    message(STATUS "GEOS_CXX_FLAGS = " ${GEOS_CXX_FLAGS})  
    message(STATUS "GEOS_EXE_LINKER_FLAGS = " ${GEOS_EXE_LINKER_FLAGS})  
    Message(STATUS "GEOS_INCLUDE_DIR = " ${GEOS_INCLUDE_DIR})
    Message(STATUS "GEOS_LINK_DIR = " ${GEOS_LINK_DIR})
ENDIF(GEOS_FOUND)

Add_subdirectory (src)
Add_subdirectory (utilities)

ENABLE_TESTING()

Set(frlen 1200,1800 CACHE STRING "Frank-Read source length" )
Set(cubel 3000 CACHE STRING "cube length" )
Set(maxSeg 200 CACHE STRING "maximum dislocation segment length" )
Set(nfrs 5 CACHE STRING "number of Frank-Read sources" )
Set(angle 1.5707963 CACHE STRING "rotation angle of TB plane" )
Set(maxDensity 1E13 CACHE STRING "maximum dislocation density" )

if(NOT seed)
String(RANDOM LENGTH 8 ALPHABET 0123456789 seed)
Message("** Random seed for test is ${seed}")
else(NOT seed)
message("** seed is defined as ${seed}")
endif(NOT seed)

Set(commands_gen_bicrystal -type 3 -frlen ${frlen} -cubel ${cubel} -maxseg ${maxSeg} -nfrsrcs ${nfrs} -seed ${seed} 
    -polycrystal 1 -nGrain 2 -pbc 7 -maxDensity ${maxDensity}
    -rotate ${angle} 
    -twinRatio 0.5 -heightRatio 2 -deltaNum 0 -nLayers 1) 

ADD_TEST(gen_bicrystal utilities/${PARADISGEN_EXE} ${commands_gen_bicrystal})

Set(commands_repart -infile paradis.data -cells 4 -domains 1,1,4 -outfile bicrystal.data)
Add_Test(repart_gen utilities/paradisrepart ${commands_repart})

Set(mpirun_commands -n 4 "./src/${PARADIS_EXE}")
if(OPENMP_MODE)
    Set(paradis_commands -n 1 -d ${PARADIS_BINARY_DIR}/bicrystal.data ${PARADIS_SOURCE_DIR}/tests/bicrystal.ctrl)
else(OPENMP_MODE)
    Set(paradis_commands -d ${PARADIS_BINARY_DIR}/bicrystal.data ${PARADIS_SOURCE_DIR}/tests/bicrystal.ctrl)
endif(OPENMP_MODE)
Add_Test(run_paradis_bicrystal "mpirun" ${mpirun_commands} ${paradis_commands})

if(0)
message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
message("CMAKE_C_FLAGS_DEBUG = ${CMAKE_C_FLAGS_DEBUG}")
message("CMAKE_C_FLAGS_MINSIZEREL = ${CMAKE_C_FLAGS_MINSIZEREL}")
message("CMAKE_C_FLAGS_RELEASE = ${CMAKE_C_FLAGS_RELEASE}")
message("CMAKE_C_FLAGS_RELWITHDEBINFO = ${CMAKE_C_FLAGS_RELWITHDEBINFO}")

message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message("CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_CXX_FLAGS_MINSIZEREL = ${CMAKE_CXX_FLAGS_MINSIZEREL}")
message("CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")
message("CMAKE_CXX_FLAGS_RELWITHDEBINFO = ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

message("CMAKE_EXE_LINKER_FLAGS = ${CMAKE_EXE_LINKER_FLAGS}")
message("CMAKE_EXE_LINKER_FLAGS_DEBUG = ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
message("CMAKE_EXE_LINKER_FLAGS_MINSIZEREL = ${CMAKE_EXE_LINKER_FLAGS_MINSIZEREL}")
message("CMAKE_EXE_LINKER_FLAGS_RELEASE = ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
message("CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO = ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}")

message("CMAKE_MODULE_LINKER_FLAGS = ${CMAKE_MODULE_LINKER_FLAGS}")
message("CMAKE_MODULE_LINKER_FLAGS_DEBUG = ${CMAKE_MODULE_LINKER_FLAGS_DEBUG}")
message("CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL = ${CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL}")
message("CMAKE_MODULE_LINKER_FLAGS_RELEASE = ${CMAKE_MODULE_LINKER_FLAGS_RELEASE}")
message("CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO = ${CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO}")

message("CMAKE_SHARED_LINKER_FLAGS = ${CMAKE_SHARED_LINKER_FLAGS}")
message("CMAKE_SHARED_LINKER_FLAGS_DEBUG = ${CMAKE_SHARED_LINKER_FLAGS_DEBUG}")
message("CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL = ${CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL}")
message("CMAKE_SHARED_LINKER_FLAGS_RELEASE = ${CMAKE_SHARED_LINKER_FLAGS_RELEASE}")
message("CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO = ${CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO}")

message("CMAKE_STATIC_LINKER_FLAGS = ${CMAKE_STATIC_LINKER_FLAGS}")
message("CMAKE_STATIC_LINKER_FLAGS_DEBUG = ${CMAKE_STATIC_LINKER_FLAGS_DEBUG}")
message("CMAKE_STATIC_LINKER_FLAGS_MINSIZEREL = ${CMAKE_STATIC_LINKER_FLAGS_MINSIZEREL}")
message("CMAKE_STATIC_LINKER_FLAGS_RELEASE = ${CMAKE_STATIC_LINKER_FLAGS_RELEASE}")
message("CMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO = ${CMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO}")

endif()
