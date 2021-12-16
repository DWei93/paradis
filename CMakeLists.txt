Project (PARADIS)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_VERBOSE_MAKEFILE OFF CACHE STRING "verbose makefile")
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

message ("** PARADIS_BINARY_DIR =  ${PARADIS_BINARY_DIR}")
message ("** PARADIS_SOURCE_DIR =  ${PARADIS_SOURCE_DIR}")


find_package(MPI REQUIRED)
if(MPI_FOUND)  
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_EXE_LINKER_FLAGS}")  
endif()  

Set(CMAKE_C_FLAGS_GDB "-ggdb -DMPI_DEBUG" CACHE STRING "C flag for gdb" )
Set(CMAKE_CXX_FLAGS_GDB "-ggdb -DMPI_DEBUG" CACHE STRING "CXX flag for gdb" )
if(NOT CMAKE_BUILD_TYPE)
    Set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE)
Set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "GDB")
message("** CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

if(${CMAKE_BUILD_TYPE} EQUAL "GDB")
    Set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_GDB} ")
    Set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_CXX_FLAGS_GDB} ")
endif()
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
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")  
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")  
        message("** Find OpenMP: OpenMP_C_FLAGS = ${OpenMP_C_FLAGS}; OpenMP_CXX_FLAGS = ${OpenMP_CXX_FLAGS}; OpenMP_EXE_LINKER_FLAGS = ${OpenMP_EXE_LINKER_FLAGS}")
    else(OPENMP_FOUND)
        Option(OPENMP_MODE OFF)
    endif(OPENMP_FOUND)  
endif(OPENMP_MODE)

if(0)
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
endif()

Add_subdirectory (src)
Add_subdirectory (utilities)

ENABLE_TESTING()

Set(frlen 1500,3000 CACHE STRING "Frank-Read source length" )
Set(cubel 4000 CACHE STRING "cube length" )
Set(maxSeg 195 CACHE STRING "maximum dislocation segment length" )
Set(nfrs 5 CACHE STRING "number of Frank-Read sources" )
Set(angle 1.5707963 CACHE STRING "rotation angle of TB plane" )
Set(maxDensity 1E13 CACHE STRING "maximum dislocation density" )

if(NOT seed)
String(RANDOM LENGTH 8 ALPHABET 0123456789 seed)
Message("** Random seed for test is ${seed}")
else(NOT seed)
message("** seed is defined as ${seed}")
endif(NOT seed)

Set( commands_gen_bicrystal
    -type 3 -frlen ${frlen} -pbc 7 -cubel ${cubel} -maxseg ${maxSeg} -nfrsrcs ${nfrs}
    -polycrystal 1 -nGrain 2 -pbc 7 -maxDensity ${maxDensity} -rotate ${angle} 
    -twinRatio 0.5 -heightRatio 2 -deltaNum 0 -nLayers 1 -nTDs 1 -seed ${seed}
    -xMatDir -1,1,0 -yMatDir -1,-1,2 -srcLenType 0 )
Set(commands_repart_bicrystal -infile paradis.data -cells 4 -domains 1,1,4 -outfile bicrystal.data)
if(NOT OPENMP_MODE OR NOT OPENMP_FOUND)
    Set(paradis_commands_bicrystal -n 1 -d ${PARADIS_BINARY_DIR}/bicrystal.data ${PARADIS_SOURCE_DIR}/tests/bicrystal.ctrl)
else(NOT OPENMP_MODE OR NOT OPENMP_FOUND)
    Set(paradis_commands_bicrystal -n 2 -d ${PARADIS_BINARY_DIR}/bicrystal.data ${PARADIS_SOURCE_DIR}/tests/bicrystal.ctrl)
endif(NOT OPENMP_MODE OR NOT OPENMP_FOUND)

Add_test(test_bicrystal-1
     utilities/${PARADISGEN_EXE} ${commands_gen_bicrystal})
Add_test(test_bicrystal-2
      utilities/paradisrepart ${commands_repart_bicrystal} )
Add_test(test_bicrystal-3
      mpirun -n 4 ./src/${PARADIS_EXE} ${paradis_commands_bicrystal})

Set(commands_gen_multilayer
    -type 3 -frlen 3 -cubel ${cubel} -maxseg ${maxSeg} -nfrsrcs ${nfrs}  -tdDensity -1  -polycrystal 1  -seed ${seed}
    -nGrain 4 -pbc 0 -nTDs 1 -maxDensity ${maxDensity} -gb 1 -loadDir 90.0,90.0  -tdLength 978 -kinkSize 1,5 )
Set(commands_repart_multilayer -infile paradis.data -cells 4 -domains 1,1,4 -outfile multilayer.data)
if(NOT OPENMP_MODE OR NOT OPENMP_FOUND)
    Set(paradis_commands_multilayer -n 1 -f ${PARADIS_SOURCE_DIR}/input/fm-ctab.Cu.m2.t5.dat -d ${PARADIS_BINARY_DIR}/multilayer.data ${PARADIS_SOURCE_DIR}/tests/multilayer.ctrl)
else(NOT OPENMP_MODE OR NOT OPENMP_FOUND)
    Set(paradis_commands_multilayer -n 2 -f ${PARADIS_SOURCE_DIR}/input/fm-ctab.Cu.m2.t5.dat -d ${PARADIS_BINARY_DIR}/multilayer.data ${PARADIS_SOURCE_DIR}/tests/multilayer.ctrl)
endif(NOT OPENMP_MODE OR NOT OPENMP_FOUND)

Add_test(test_multilayer-1
     utilities/${PARADISGEN_EXE} ${commands_gen_multilayer})
Add_test(test_multilayer-2
      utilities/paradisrepart ${commands_repart_multilayer} )
Add_test(test_multilayer-3
      mpirun -n 4 ./src/${PARADIS_EXE} ${paradis_commands_multilayer})

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
