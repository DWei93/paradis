Set ( PARADISGEN_SRCS 
               ParadisGen.cpp   
               CreateConfig.cpp      
               InitRemesh.cpp
               ../src/Decomp.cpp            
               ../src/DLBfreeOld.cpp        
               ../src/FindPreciseGlidePlane.cpp  
               ../src/FMComm.cpp           
               ../src/FMSigma2.cpp          
               ../src/FMSupport.cpp         
               ../src/FreeInitArrays.cpp    
               ../src/GetNewNativeNode.cpp  
               ../src/Heap.cpp              
               ../src/InitCellDomains.cpp   
               ../src/InitCellNatives.cpp   
               ../src/InitCellNeighbors.cpp 
               ../src/InitSendDomains.cpp   
               ../src/InitRemoteDomains.cpp 
               ../src/Matrix.cpp            
               ../src/MemCheck.cpp          
               ../src/Meminfo.cpp         
               ../src/Param.cpp             
               ../src/Parse.cpp             
               ../src/PickScrewGlidePlane.cpp 
               ../src/QueueOps.cpp          
               ../src/RBDecomp.cpp          
               ../src/ReadRestart.cpp       
               ../src/ReadBinaryRestart.cpp 
               ../src/RSDecomp.cpp          
               ../src/Timer.cpp             
               ../src/Util.cpp              
               ../src/Math.cpp              
               ../src/Grain.cpp )

Set ( PARADISREPART_SRCS
                    ParadisRepart.cpp  
                    ../src/Decomp.cpp            
                    ../src/DLBfreeOld.cpp        
                    ../src/FindPreciseGlidePlane.cpp  
                    ../src/FMComm.cpp            
                    ../src/FMSigma2.cpp          
                    ../src/FMSupport.cpp         
                    ../src/FreeInitArrays.cpp    
                    ../src/Heap.cpp              
                    ../src/InitCellDomains.cpp   
                    ../src/InitCellNatives.cpp   
                    ../src/InitCellNeighbors.cpp 
                    ../src/InitRemoteDomains.cpp 
                    ../src/InitSendDomains.cpp   
                    ../src/Matrix.cpp            
                    ../src/MemCheck.cpp          
                    ../src/Param.cpp             
                    ../src/Parse.cpp             
                    ../src/PickScrewGlidePlane.cpp 
                    ../src/QueueOps.cpp          
                    ../src/ReadRestart.cpp       
                    ../src/ReadBinaryRestart.cpp 
                    ../src/RBDecomp.cpp          
                    ../src/RSDecomp.cpp          
                    ../src/Timer.cpp             
                    ../src/Util.cpp              
                    ../src/Math.cpp              
                    ../src/WriteRestart.cpp		
					../src/Grain.cpp )
Set( CALCDENSITY_SRCS
                    CalcDensity.cpp      
                    ../src/Decomp.cpp            
                    ../src/DLBfreeOld.cpp        
                    ../src/FindPreciseGlidePlane.cpp  
                    ../src/FMComm.cpp            
                    ../src/FMSigma2.cpp          
                    ../src/FMSupport.cpp         
                    ../src/FreeInitArrays.cpp    
                    ../src/Heap.cpp              
                    ../src/InitCellDomains.cpp   
                    ../src/InitCellNatives.cpp   
                    ../src/InitCellNeighbors.cpp 
                    ../src/InitRemoteDomains.cpp 
                    ../src/InitSendDomains.cpp   
                    ../src/Matrix.cpp            
                    ../src/MemCheck.cpp         
                    ../src/Param.cpp             
                    ../src/Parse.cpp             
                    ../src/PickScrewGlidePlane.cpp 
                    ../src/QueueOps.cpp          
                    ../src/ReadRestart.cpp       
                    ../src/ReadBinaryRestart.cpp 
                    ../src/RBDecomp.cpp          
                    ../src/RSDecomp.cpp          
                    ../src/Timer.cpp             
                    ../src/Util.cpp              
                    ../src/Math.cpp              
					../src/Grain.cpp )

Set ( STRESSTABLEGEN_SRCS
                      StressTableGen.cpp 
                      ../src/FindPreciseGlidePlane.cpp 
                      ../src/Heap.cpp           
                      ../src/Matrix.cpp         
                      ../src/MemCheck.cpp       
                      ../src/PickScrewGlidePlane.cpp 
                      ../src/QueueOps.cpp       
                      ../src/Util.cpp      		
                      ../src/Math.cpp      		
					  ../src/Grain.cpp )

Set( PARADISCONVERT_SRCS 
                     ParadisConvert.cpp 
                     ../src/Decomp.cpp            
                     ../src/DLBfreeOld.cpp        
                     ../src/FindPreciseGlidePlane.cpp  
                     ../src/FMComm.cpp            
                     ../src/FMSigma2.cpp          
                     ../src/FMSupport.cpp         
                     ../src/FreeInitArrays.cpp    
                     ../src/Heap.cpp              
                     ../src/InitCellDomains.cpp   
                     ../src/InitCellNatives.cpp  
                     ../src/InitCellNeighbors.cpp 
                     ../src/InitRemoteDomains.cpp 
					 ../src/InitSendDomains.cpp   
                     ../src/Matrix.cpp            
                     ../src/MemCheck.cpp          
                     ../src/Param.cpp             
                     ../src/Parse.cpp             
                     ../src/PickScrewGlidePlane.cpp 
                     ../src/QueueOps.cpp          
                     ../src/ReadRestart.cpp       
                     ../src/ReadBinaryRestart.cpp 
                     ../src/RBDecomp.cpp          
                     ../src/RSDecomp.cpp          
                     ../src/Timer.cpp             
                     ../src/Util.cpp         	
                     ../src/Math.cpp         	
					 ../src/Grain.cpp )

Set(CTABLEGEN_MAIN_SRC  ../src/CTableGen.cpp
                        ../src/CorrectionTable.cpp       
                        ../src/FindPreciseGlidePlane.cpp 
                        ../src/FMSigma2.cpp              
                        ../src/FMSupport.cpp             
                        ../src/Matrix.cpp                
                        ../src/Heap.cpp                  
                        ../src/MemCheck.cpp              
                        ../src/PickScrewGlidePlane.cpp   
                        ../src/QueueOps.cpp              
                        ../src/Util.cpp					 
                        ../src/Math.cpp	                 
                        ../src/Grain.cpp	 )


IF(GEOS_FOUND)
INCLUDE_DIRECTORIES(../include ${GEOS_INCLUDE_DIR})    
LINK_DIRECTORIES(${GEOS_LINK_DIR})
ELSE(GEOS_FOUND)
INCLUDE_DIRECTORIES(../include)
ENDIF(GEOS_FOUND)

add_definitions(-DNO_XWINDOW)
if(MPI_FOUND)
    add_definitions(-DDPARALLEL=1)
endif(MPI_FOUND)

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

