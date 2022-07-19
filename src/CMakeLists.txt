Set ( PARADIS_SRCS 
#bigtest.h   GeometryTestFactory.cc TestSweepLineSpeed.cc  
#bci.c   tdspl.c tools.c
CppFunctions.cc
CellCharge.cc               FindPreciseGlidePlane.cc   InputSanity.cc            ParadisThread.cc            SplitSurfaceNodes.cc
Collision.cc                FixRemesh.cc               LoadCurve.cc              Param.cc                    Tecplot.cc
CommSendGhostPlanes.cc      FMComm.cc                  LocalSegForces.cc         Parse.cc                    Timer.cc
CommSendGhosts.cc           FMSigma2.cc                Main.cc                   PickScrewGlidePlane.cc      Topology.cc
CommSendMirrorNodes.cc      FMSupport.cc                                          Plot.cc                     TrapezoidIntegrator.cc
CommSendRemesh.cc           ForwardEulerIntegrator.cc                             PredictiveCollision.cc      Util.cc
CommSendSecondaryGhosts.cc  FreeInitArrays.cc          Matrix.cc                 ProximityCollision.cc       WriteArms.cc
CommSendSegments.cc         GenerateOutput.cc          MemCheck.cc               QueueOps.cc                 WriteBinaryRestart.cc
CommSendVelocity.cc         GetDensityDelta.cc         Meminfo.cc                RBDecomp.cc                 WriteDensFlux.cc
CorrectionTable.cc          GetNewGhostNode.cc         Migrate.cc                ReadBinaryRestart.cc        WriteDensityField.cc
CrossSlipBCC.cc             GetNewNativeNode.cc        MobilityLaw_BCC_0b.cc     ReadRestart.cc              WriteForce.cc
CrossSlip.cc                Gnuplot.cc                 MobilityLaw_BCC_0.cc      RemapInitialTags.cc         WriteFragments.cc
CrossSlipFCC.cc             GrainBoundary.cc           MobilityLaw_BCC_glide.cc  Remesh.cc                   WritePoleFig.cc
                             Grain.cc                   MobilityLaw_FCC_0b.cc     RemeshRule_2.cc             WritePovray.cc
DebugFunctions.cc           Heap.cc                    MobilityLaw_FCC_0.cc      RemeshRule_3.cc             WriteProp.cc
Decomp.cc                   InitCellDomains.cc         MobilityLaw_FCC_climb.cc  RemoteSegForces.cc          WriteRestart.cc
DeltaPlasticStrain_BCC.cc   InitCellNatives.cc         NodeForce.cc              RemoveNode.cc               WriteVelocity.cc
DeltaPlasticStrain.cc       InitCellNeighbors.cc       NodeVelocity.cc           ResetGlidePlanes.cc         WriteVisit.cc
DeltaPlasticStrain_FCC.cc   InitHome.cc                OsmoticForce.cc           RSDecomp.cc
deWitInteraction.cc         Initialize.cc              ParadisFinish.cc          SemiInfiniteSegSegForce.cc
DisableUnneededParams.cc    InitRemoteDomains.cc       ParadisInit.cc            SortNativeNodes.cc
DLBfreeOld.cc               InitSendDomains.cc         ParadisStep.cc            SortNodesForCollision.cc    Math.cc 
MobilityLaw_FCC_Nolinear.cc MobilityLaw_Wei.cc )
#      DisplayC.cc display.cc )

Set (CTABLEGEN_SRCS 
                 CTableGen.cc 
                 CorrectionTable.cc       
                 FindPreciseGlidePlane.cc 
                 FMSigma2.cc              
                 FMSupport.cc             
                 Matrix.cc                
                 Heap.cc                  
                 MemCheck.cc              
                 PickScrewGlidePlane.cc   
                 QueueOps.cc              
                 Util.cc					
                 Math.cc					
		         Grain.cc)

if(GEOS_FOUND)
INCLUDE_DIRECTORIES(../include ${MPI_INCLUDE_PATH} ${GEOS_INCLUDE_DIR} )
else(GEOS_FOUND)
INCLUDE_DIRECTORIES(../include ${MPI_INCLUDE_PATH})
endif(GEOS_FOUND)

if(MPI_FOUND)
add_definitions(-DNO_XWINDOW -DPARALLEL=1)
else(MPI_FOUND)
add_definitions(-DNO_XWINDOW)
endif(MPI_FOUND)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas -fopenmp" )
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp " )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")

Option(JOG_PINNING "Build ${PARADIS_EXE} using definition JOG_PINNING" OFF)
if(JOG_PINNING)
    add_definitions(-DJOG_PINNING)
endif()
Message("** Macro JOG_PINNING is ${JOG_PINNING}")

Option(VELOCITY_STATISTICS "Build ${PARADIS_EXE} using definition VELOCITY_STATISTICS" OFF)
if(VELOCITY_STATISTICS)
    add_definitions(-DVELOCITY_STATISTICS)
endif()
Message("** Macro VELOCITY_STATISTICS is ${VELOCITY_STATISTICS}")

Option(NO_GB_BARRIER "Build ${PARADIS_EXE} using definition NO_GB_BARRIER" OFF)
if(NO_GB_BARRIER)
    Message(STATUS "Macro NO_GB_BARRIER is on")
    add_definitions(-DNO_GB_BARRIER)
endif(NO_GB_BARRIER)

Option(DEBUG_PLASTIC_STRAIN "debug plastic strain field" OFF)
if(DEBUG_PLASTIC_STRAIN)
    add_definitions(-DDEBUG_PLASTIC_STRAIN)
endif()

Add_executable(${PARADIS_EXE} ${PARADIS_SRCS})
target_link_libraries(${PARADIS_EXE} ${MPI_CXX_LIBRARIES}) # this must be  MPI_CXX_LIBRARIES for cpp file
Add_definitions(-DSIMLPLE_REACTION_ON_GB)


Install(TARGETS ${PARADIS_EXE} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
