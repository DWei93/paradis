Set ( PARADISGEN_SRCS 
               ParadisGen.C   
               CreateConfig.C      
               InitRemesh.C
               ../src/Decomp.C            
               ../src/DLBfreeOld.C        
               ../src/FindPreciseGlidePlane.C  
               ../src/FMComm.C           
               ../src/FMSigma2.C          
               ../src/FMSupport.C         
               ../src/FreeInitArrays.C    
               ../src/GetNewNativeNode.C  
               ../src/Heap.C              
               ../src/InitCellDomains.C   
               ../src/InitCellNatives.C   
               ../src/InitCellNeighbors.C 
               ../src/InitSendDomains.C   
               ../src/InitRemoteDomains.C 
               ../src/Matrix.C            
               ../src/MemCheck.C          
               ../src/Meminfo.C         
               ../src/Param.C             
               ../src/Parse.C             
               ../src/PickScrewGlidePlane.C 
               ../src/QueueOps.C          
               ../src/RBDecomp.C          
               ../src/ReadRestart.C       
               ../src/ReadBinaryRestart.C 
               ../src/RSDecomp.C          
               ../src/Timer.C             
               ../src/Util.C              
               ../src/Math.C              
               ../src/Grain.C )

Set ( PARADISREPART_SRCS
                    ParadisRepart.C  
                    ../src/Decomp.C            
                    ../src/DLBfreeOld.C        
                    ../src/FindPreciseGlidePlane.C  
                    ../src/FMComm.C            
                    ../src/FMSigma2.C          
                    ../src/FMSupport.C         
                    ../src/FreeInitArrays.C    
                    ../src/Heap.C              
                    ../src/InitCellDomains.C   
                    ../src/InitCellNatives.C   
                    ../src/InitCellNeighbors.C 
                    ../src/InitRemoteDomains.C 
                    ../src/InitSendDomains.C   
                    ../src/Matrix.C            
                    ../src/MemCheck.C          
                    ../src/Param.C             
                    ../src/Parse.C             
                    ../src/PickScrewGlidePlane.C 
                    ../src/QueueOps.C          
                    ../src/ReadRestart.C       
                    ../src/ReadBinaryRestart.C 
                    ../src/RBDecomp.C          
                    ../src/RSDecomp.C          
                    ../src/Timer.C             
                    ../src/Util.C              
                    ../src/Math.C              
                    ../src/WriteRestart.C		
					../src/Grain.C )
Set( CALCDENSITY_SRCS
                    CalcDensity.C      
                    ../src/Decomp.C            
                    ../src/DLBfreeOld.C        
                    ../src/FindPreciseGlidePlane.C  
                    ../src/FMComm.C            
                    ../src/FMSigma2.C          
                    ../src/FMSupport.C         
                    ../src/FreeInitArrays.C    
                    ../src/Heap.C              
                    ../src/InitCellDomains.C   
                    ../src/InitCellNatives.C   
                    ../src/InitCellNeighbors.C 
                    ../src/InitRemoteDomains.C 
                    ../src/InitSendDomains.C   
                    ../src/Matrix.C            
                    ../src/MemCheck.C         
                    ../src/Param.C             
                    ../src/Parse.C             
                    ../src/PickScrewGlidePlane.C 
                    ../src/QueueOps.C          
                    ../src/ReadRestart.C       
                    ../src/ReadBinaryRestart.C 
                    ../src/RBDecomp.C          
                    ../src/RSDecomp.C          
                    ../src/Timer.C             
                    ../src/Util.C              
                    ../src/Math.C              
					../src/Grain.C )

Set ( STRESSTABLEGEN_SRCS
                      StressTableGen.C 
                      ../src/FindPreciseGlidePlane.C 
                      ../src/Heap.C           
                      ../src/Matrix.C         
                      ../src/MemCheck.C       
                      ../src/PickScrewGlidePlane.C 
                      ../src/QueueOps.C       
                      ../src/Util.C      		
                      ../src/Math.C      		
					  ../src/Grain.C )

Set( PARADISCONVERT_SRCS 
                     ParadisConvert.C 
                     ../src/Decomp.C            
                     ../src/DLBfreeOld.C        
                     ../src/FindPreciseGlidePlane.C  
                     ../src/FMComm.C            
                     ../src/FMSigma2.C          
                     ../src/FMSupport.C         
                     ../src/FreeInitArrays.C    
                     ../src/Heap.C              
                     ../src/InitCellDomains.C   
                     ../src/InitCellNatives.C  
                     ../src/InitCellNeighbors.C 
                     ../src/InitRemoteDomains.C 
					 ../src/InitSendDomains.C   
                     ../src/Matrix.C            
                     ../src/MemCheck.C          
                     ../src/Param.C             
                     ../src/Parse.C             
                     ../src/PickScrewGlidePlane.C 
                     ../src/QueueOps.C          
                     ../src/ReadRestart.C       
                     ../src/ReadBinaryRestart.C 
                     ../src/RBDecomp.C          
                     ../src/RSDecomp.C          
                     ../src/Timer.C             
                     ../src/Util.C         	
                     ../src/Math.C         	
					 ../src/Grain.C )

Set(CTABLEGEN_MAIN_SRC  ../src/CTableGen.C
                        ../src/CorrectionTable.C       
                        ../src/FindPreciseGlidePlane.C 
                        ../src/FMSigma2.C              
                        ../src/FMSupport.C             
                        ../src/Matrix.C                
                        ../src/Heap.C                  
                        ../src/MemCheck.C              
                        ../src/PickScrewGlidePlane.C   
                        ../src/QueueOps.C              
                        ../src/Util.C					 
                        ../src/Math.C	                 
                        ../src/Grain.C	 )


IF(GEOS_FOUND)
INCLUDE_DIRECTORIES(../include ${GEOS_INCLUDE_DIR})    
LINK_DIRECTORIES(${GEOS_LINK_DIR})
ELSE(GEOS_FOUND)
INCLUDE_DIRECTORIES(../include)
ENDIF(GEOS_FOUND)

add_definitions(-DNO_XWINDOW )

add_executable(${PARADISGEN_EXE} ${PARADISGEN_SRCS}     )
add_executable(paradisrepart     ${PARADISREPART_SRCS}  )
add_executable(calcdensity       ${CALCDENSITY_SRCS}    )
add_executable(stresstablegen    ${STRESSTABLEGEN_SRCS} )
add_executable(paradisconvert    ${PARADISCONVERT_SRCS} )
add_executable(ctablegen         ${CTABLEGEN_MAIN_SRC}  )

target_link_libraries(${PARADISGEN_EXE} ${MPI_CXX_LIBRARIES}) 
target_link_libraries(paradisrepart     ${MPI_CXX_LIBRARIES})
target_link_libraries(calcdensity       ${MPI_CXX_LIBRARIES}) 
target_link_libraries(stresstablegen    ${MPI_CXX_LIBRARIES}) 
target_link_libraries(paradisconvert    ${MPI_CXX_LIBRARIES})
target_link_libraries(ctablegen         ${MPI_CXX_LIBRARIES})

Install(TARGETS ${PARADISGEN_EXE} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
Install(TARGETS paradisrepart     RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
Install(TARGETS calcdensity       RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
Install(TARGETS stresstablegen    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
Install(TARGETS paradisconvert    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
Install(TARGETS ctablegen         RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

