Set ( PARADISGEN_SRCS 
               ParadisGen.cc   
               CreateConfig.cc      
               InitRemesh.cc
               ../src/Decomp.cc            
               ../src/DLBfreeOld.cc        
               ../src/FindPreciseGlidePlane.cc  
               ../src/FMComm.cc           
               ../src/FMSigma2.cc          
               ../src/FMSupport.cc         
               ../src/FreeInitArrays.cc    
               ../src/GetNewNativeNode.cc  
               ../src/Heap.cc              
               ../src/InitCellDomains.cc   
               ../src/InitCellNatives.cc   
               ../src/InitCellNeighbors.cc 
               ../src/InitSendDomains.cc   
               ../src/InitRemoteDomains.cc 
               ../src/Matrix.cc            
               ../src/MemCheck.cc          
               ../src/Meminfo.cc         
               ../src/Param.cc             
               ../src/Parse.cc             
               ../src/PickScrewGlidePlane.cc 
               ../src/QueueOps.cc          
               ../src/RBDecomp.cc          
               ../src/ReadRestart.cc       
               ../src/ReadBinaryRestart.cc 
               ../src/RSDecomp.cc          
               ../src/Timer.cc             
               ../src/Util.cc              
               ../src/Math.cc              
               ../src/Grain.cc )

Set ( PARADISREPART_SRCS
                    ParadisRepart.cc  
                    ../src/Decomp.cc            
                    ../src/DLBfreeOld.cc        
                    ../src/FindPreciseGlidePlane.cc  
                    ../src/FMComm.cc            
                    ../src/FMSigma2.cc          
                    ../src/FMSupport.cc         
                    ../src/FreeInitArrays.cc    
                    ../src/Heap.cc              
                    ../src/InitCellDomains.cc   
                    ../src/InitCellNatives.cc   
                    ../src/InitCellNeighbors.cc 
                    ../src/InitRemoteDomains.cc 
                    ../src/InitSendDomains.cc   
                    ../src/Matrix.cc            
                    ../src/MemCheck.cc          
                    ../src/Param.cc             
                    ../src/Parse.cc             
                    ../src/PickScrewGlidePlane.cc 
                    ../src/QueueOps.cc          
                    ../src/ReadRestart.cc       
                    ../src/ReadBinaryRestart.cc 
                    ../src/RBDecomp.cc          
                    ../src/RSDecomp.cc          
                    ../src/Timer.cc             
                    ../src/Util.cc              
                    ../src/Math.cc              
                    ../src/WriteRestart.cc		
					../src/Grain.cc )
Set( CALCDENSITY_SRCS
                    CalcDensity.cc      
                    ../src/Decomp.cc            
                    ../src/DLBfreeOld.cc        
                    ../src/FindPreciseGlidePlane.cc  
                    ../src/FMComm.cc            
                    ../src/FMSigma2.cc          
                    ../src/FMSupport.cc         
                    ../src/FreeInitArrays.cc    
                    ../src/Heap.cc              
                    ../src/InitCellDomains.cc   
                    ../src/InitCellNatives.cc   
                    ../src/InitCellNeighbors.cc 
                    ../src/InitRemoteDomains.cc 
                    ../src/InitSendDomains.cc   
                    ../src/Matrix.cc            
                    ../src/MemCheck.cc         
                    ../src/Param.cc             
                    ../src/Parse.cc             
                    ../src/PickScrewGlidePlane.cc 
                    ../src/QueueOps.cc          
                    ../src/ReadRestart.cc       
                    ../src/ReadBinaryRestart.cc 
                    ../src/RBDecomp.cc          
                    ../src/RSDecomp.cc          
                    ../src/Timer.cc             
                    ../src/Util.cc              
                    ../src/Math.cc              
					../src/Grain.cc )

Set ( STRESSTABLEGEN_SRCS
                      StressTableGen.cc 
                      ../src/FindPreciseGlidePlane.cc 
                      ../src/Heap.cc           
                      ../src/Matrix.cc         
                      ../src/MemCheck.cc       
                      ../src/PickScrewGlidePlane.cc 
                      ../src/QueueOps.cc       
                      ../src/Util.cc      		
                      ../src/Math.cc      		
					  ../src/Grain.cc )

Set( PARADISCONVERT_SRCS 
                     ParadisConvert.cc 
                     ../src/Decomp.cc            
                     ../src/DLBfreeOld.cc        
                     ../src/FindPreciseGlidePlane.cc  
                     ../src/FMComm.cc            
                     ../src/FMSigma2.cc          
                     ../src/FMSupport.cc         
                     ../src/FreeInitArrays.cc    
                     ../src/Heap.cc              
                     ../src/InitCellDomains.cc   
                     ../src/InitCellNatives.cc  
                     ../src/InitCellNeighbors.cc 
                     ../src/InitRemoteDomains.cc 
					 ../src/InitSendDomains.cc   
                     ../src/Matrix.cc            
                     ../src/MemCheck.cc          
                     ../src/Param.cc             
                     ../src/Parse.cc             
                     ../src/PickScrewGlidePlane.cc 
                     ../src/QueueOps.cc          
                     ../src/ReadRestart.cc       
                     ../src/ReadBinaryRestart.cc 
                     ../src/RBDecomp.cc          
                     ../src/RSDecomp.cc          
                     ../src/Timer.cc             
                     ../src/Util.cc         	
                     ../src/Math.cc         	
					 ../src/Grain.cc )

Set(CTABLEGEN_MAIN_SRC  ../src/CTableGen.cc
                        ../src/CorrectionTable.cc       
                        ../src/FindPreciseGlidePlane.cc 
                        ../src/FMSigma2.cc              
                        ../src/FMSupport.cc             
                        ../src/Matrix.cc                
                        ../src/Heap.cc                  
                        ../src/MemCheck.cc              
                        ../src/PickScrewGlidePlane.cc   
                        ../src/QueueOps.cc              
                        ../src/Util.cc					 
                        ../src/Math.cc	                 
                        ../src/Grain.cc	 )


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

