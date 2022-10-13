Set (PARADIS_SRCS 
CppFunctions.cpp             
CellCharge.cpp               FindPreciseGlidePlane.cpp   InputSanity.cpp            ParadisThread.cpp            SplitSurfaceNodes.cpp
Collision.cpp                FixRemesh.cpp               LoadCurve.cpp              Param.cpp                    Tecplot.cpp
CommSendGhostPlanes.cpp      FMComm.cpp                  LocalSegForces.cpp         Parse.cpp                    Timer.cpp
CommSendGhosts.cpp           FMSigma2.cpp                Main.cpp                   PickScrewGlidePlane.cpp      Topology.cpp
CommSendMirrorNodes.cpp      FMSupport.cpp                                          Plot.cpp                     TrapezoidIntegrator.cpp
CommSendRemesh.cpp           ForwardEulerIntegrator.cpp                             PredictiveCollision.cpp      Util.cpp
CommSendSecondaryGhosts.cpp  FreeInitArrays.cpp          Matrix.cpp                 ProximityCollision.cpp       WriteArms.cpp
CommSendSegments.cpp         GenerateOutput.cpp          MemCheck.cpp               QueueOps.cpp                 WriteBinaryRestart.cpp
CommSendVelocity.cpp         GetDensityDelta.cpp         Meminfo.cpp                RBDecomp.cpp                 WriteDensFlux.cpp
CorrectionTable.cpp          GetNewGhostNode.cpp         Migrate.cpp                ReadBinaryRestart.cpp        WriteDensityField.cpp
CrossSlipBCC.cpp             GetNewNativeNode.cpp        MobilityLaw_BCC_0b.cpp     ReadRestart.cpp              WriteForce.cpp
CrossSlip.cpp                Gnuplot.cpp                 MobilityLaw_BCC_0.cpp      RemapInitialTags.cpp         WriteFragments.cpp
CrossSlipFCC.cpp             GrainBoundary.cpp           MobilityLaw_BCC_glide.cpp  Remesh.cpp                   WritePoleFig.cpp
                             Grain.cpp                   MobilityLaw_FCC_0b.cpp     RemeshRule_2.cpp             WritePovray.cpp
DebugFunctions.cpp           Heap.cpp                    MobilityLaw_FCC_0.cpp      RemeshRule_3.cpp             WriteProp.cpp
Decomp.cpp                   InitCellDomains.cpp         MobilityLaw_FCC_climb.cpp  RemoteSegForces.cpp          WriteRestart.cpp
DeltaPlasticStrain_BCC.cpp   InitCellNatives.cpp         NodeForce.cpp              RemoveNode.cpp               WriteVelocity.cpp
DeltaPlasticStrain.cpp       InitCellNeighbors.cpp       NodeVelocity.cpp           ResetGlidePlanes.cpp         WriteVisit.cpp
DeltaPlasticStrain_FCC.cpp   InitHome.cpp                OsmoticForce.cpp           RSDecomp.cpp
deWitInteraction.cpp         Initialize.cpp              ParadisFinish.cpp          SemiInfiniteSegSegForce.cpp
DisableUnneededParams.cpp    InitRemoteDomains.cpp       ParadisInit.cpp            SortNativeNodes.cpp
DLBfreeOld.cpp               InitSendDomains.cpp         ParadisStep.cpp            SortNodesForCollision.cpp    Math.cpp 
MobilityLaw_FCC_Nolinear.cpp )


Set (CTABLEGEN_SRCS 
                 CTableGen.cpp 
                 CorrectionTable.cpp       
                 FindPreciseGlidePlane.cpp 
                 FMSigma2.cpp              
                 FMSupport.cpp             
                 Matrix.cpp                
                 Heap.cpp                  
                 MemCheck.cpp              
                 PickScrewGlidePlane.cpp   
                 QueueOps.cpp              
                 Util.cpp					
                 Math.cpp					
		         Grain.cpp)

AddDef(-DNO_XWINDOW)
if(MPI_FOUND)
    AddDef(-DPARALLEL=1)
endif(MPI_FOUND)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas -fopenmp" )
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp " )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")

Option(JOG_PINNING "Build ${PARADIS_EXE} using definition JOG_PINNING" OFF)
if(JOG_PINNING)
    AddDef(-DJOG_PINNING)
endif()

Option(VELOCITY_STATISTICS "Build ${PARADIS_EXE} using definition VELOCITY_STATISTICS" OFF)
if(VELOCITY_STATISTICS)
    AddDef(-DVELOCITY_STATISTICS)
endif()

if(PETSC_DEFINITION)
    set(PARADIS_SRCS PetscUtil.cpp PetscVelocityKernel.cpp ${PARADIS_SRCS})
    AddDef(-DPETSC)
    if(PETSC_DEBUG)
        AddDef(-DPETSC_DEBUG)
    endif()
endif()

Option(NO_GB_BARRIER "Build ${PARADIS_EXE} using definition NO_GB_BARRIER" OFF)
if(NO_GB_BARRIER)
    AddDef(-DNO_GB_BARRIER)
endif(NO_GB_BARRIER)

Option(CHECK_ALL "Breakpoint detection definition" OFF)
if(CHECK_ALL)
    AddDef(-DCHECK_ALL)
endif()

Option(DEBUG_PLASTIC_STRAIN "Debug plastic strain field" OFF)
if(DEBUG_PLASTIC_STRAIN)
    AddDef(-DDEBUG_PLASTIC_STRAIN)
endif()

Option(SIMPLE_REACTION_ON_GB "Simple GB reaction (default on)" ON)
if(SIMPLE_REACTION_ON_GB)
    AddDef(-DSIMPLE_REACTION_ON_GB)
endif()

Option(DEBUG_STEP "Debug paradis step (default off)" OFF)
if(DEBUG_STEP)
    AddDef(-DDEBUG_STEP)
endif()

add_definitions(${compile_defs})
message(STATUS "${PARADIS_EXE} compile definitions are (${compile_defs})")

Add_executable(${PARADIS_EXE} ${PARADIS_SRCS})

if(GEOS_FOUND)
target_include_directories(${PARADIS_EXE} PUBLIC  ${PETSC_INCLUDES} ../include)
else()
target_include_directories(${PARADIS_EXE} PUBLIC ${PETSC_INCLUDES} ../include ${MPI_INCLUDE_PATH})
endif()

target_link_libraries(${PARADIS_EXE} ${MPI_CXX_LIBRARIES} ${PETSC_WITH_EXTERNAL_LIB}) # this must be  MPI_CXX_LIBRARIES for cpp file


Install(TARGETS ${PARADIS_EXE} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)




