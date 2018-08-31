Set ( PARADIS_SRCS 
#bigtest.h   GeometryTestFactory.cpp TestSweepLineSpeed.cpp  
bci.c   tdspl.c tools.c
CellCharge.C               FindPreciseGlidePlane.C   InputSanity.C            ParadisThread.C            SplitSurfaceNodes.C
Collision.C                FixRemesh.C               LoadCurve.C              Param.C                    Tecplot.C
CommSendGhostPlanes.C      FMComm.C                  LocalSegForces.C         Parse.C                    Timer.C
CommSendGhosts.C           FMSigma2.C                Main.C                   PickScrewGlidePlane.C      Topology.C
CommSendMirrorNodes.C      FMSupport.C                                        Plot.C                     TrapezoidIntegrator.C
CommSendRemesh.C           ForwardEulerIntegrator.C                           PredictiveCollision.C      Util.C
CommSendSecondaryGhosts.C  FreeInitArrays.C          Matrix.C                 ProximityCollision.C       WriteArms.C
CommSendSegments.C         GenerateOutput.C          MemCheck.C               QueueOps.C                 WriteBinaryRestart.C
CommSendVelocity.C         GetDensityDelta.C         Meminfo.C                RBDecomp.C                 WriteDensFlux.C
CorrectionTable.C          GetNewGhostNode.C         Migrate.C                ReadBinaryRestart.C        WriteDensityField.C
CrossSlipBCC.C             GetNewNativeNode.C        MobilityLaw_BCC_0b.C     ReadRestart.C              WriteForce.C
CrossSlip.C                Gnuplot.C                 MobilityLaw_BCC_0.C      RemapInitialTags.C         WriteFragments.C
CrossSlipFCC.C             GrainBoundary.C           MobilityLaw_BCC_glide.C  Remesh.C                   WritePoleFig.C
                           Grain.C                   MobilityLaw_FCC_0b.C     RemeshRule_2.C             WritePovray.C
DebugFunctions.C           Heap.C                    MobilityLaw_FCC_0.C      RemeshRule_3.C             WriteProp.C
Decomp.C                   InitCellDomains.C         MobilityLaw_FCC_climb.C  RemoteSegForces.C          WriteRestart.C
DeltaPlasticStrain_BCC.C   InitCellNatives.C         NodeForce.C              RemoveNode.C               WriteVelocity.C
DeltaPlasticStrain.C       InitCellNeighbors.C       NodeVelocity.C           ResetGlidePlanes.C         WriteVisit.C
DeltaPlasticStrain_FCC.C   InitHome.C                OsmoticForce.C           RSDecomp.C
deWitInteraction.C         Initialize.C              ParadisFinish.C          SemiInfiniteSegSegForce.C
DisableUnneededParams.C    InitRemoteDomains.C       ParadisInit.C            SortNativeNodes.C
DLBfreeOld.C               InitSendDomains.C         ParadisStep.C            SortNodesForCollision.C   Math.C)
#      DisplayC.C display.C )

Set (CTABLEGEN_SRCS 
bci.c   tdspl.c tools.c
                 CTableGen.C 
                 CorrectionTable.C       
                 FindPreciseGlidePlane.C 
                 FMSigma2.C              
                 FMSupport.C             
                 Matrix.C                
                 Heap.C                  
                 MemCheck.C              
                 PickScrewGlidePlane.C   
                 QueueOps.C              
                 Util.C					
                 Math.C					
				 Grain.C)

INCLUDE_DIRECTORIES(../include ${GEOS_INCLUDE_DIR})    
LINK_DIRECTORIES(${GEOS_LINK_DIR})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-unknown-pragmas -DPARALLEL=1 " )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 " )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DNO_XPM -DNO_GENERAL -D_SEM_SEMUN_UNDEFINED")
add_executable(paradis ${PARADIS_SRCS})
add_definitions(-DNO_XWINDOW -DPARALLEL=1)


target_link_libraries(paradis m geos pthread )
Install(TARGETS paradis
        RUNTIME DESTINATION ${PARADIS_SOURCE_DIR}/bin)

