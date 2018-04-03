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
                    ../src/WriteRestart.C		
					../src/Grain.C )
Set( CALCDENSITY_SRCS
                    CalcDensity.C      
                    Decomp.C            
                    DLBfreeOld.C        
                    FindPreciseGlidePlane.C  
                    FMComm.C            
                    FMSigma2.C          
                    FMSupport.C         
                    FreeInitArrays.C    
                    Heap.C              
                    InitCellDomains.C   
                    InitCellNatives.C   
                    InitCellNeighbors.C 
                    InitRemoteDomains.C 
                    InitSendDomains.C   
                    Matrix.C            
                    MemCheck.C         
                    Param.C             
                    Parse.C             
                    PickScrewGlidePlane.C 
                    QueueOps.C          
                    ReadRestart.C       
                    ReadBinaryRestart.C 
                    RBDecomp.C          
                    RSDecomp.C          
                    Timer.C             
                    Util.C              
					Grain.C )

Set ( STRESSTABLEGEN_SRCS
                      StressTableGen.C 
                      FindPreciseGlidePlane.C 
                      Heap.C           
                      Matrix.C         
                      MemCheck.C       
                      PickScrewGlidePlane.C 
                      QueueOps.C       
                      Util.C      		
					  Grain.C )

Set( PARADISCONVERT_SRCS 
                     ParadisConvert.C 
                     Decomp.C            
                     DLBfreeOld.C        
                     FindPreciseGlidePlane.C  
                     FMComm.C            
                     FMSigma2.C          
                     FMSupport.C         
                     FreeInitArrays.C    
                     Heap.C              
                     InitCellDomains.C   
                     InitCellNatives.C  
                     InitCellNeighbors.C 
                     InitRemoteDomains.C 
					 InitSendDomains.C   
                     Matrix.C            
                     MemCheck.C          
                     Param.C             
                     Parse.C             
                     PickScrewGlidePlane.C 
                     QueueOps.C          
                     ReadRestart.C       
                     ReadBinaryRestart.C 
                     RBDecomp.C          
                     RSDecomp.C          
                     Timer.C             
                     Util.C         	
					 Grain.C )
INCLUDE_DIRECTORIES(../include ${GEOS_INCLUDE_DIR})    
LINK_DIRECTORIES(${GEOS_LINK_DIR})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-unknown-pragmas -DPARALLEL=1 " )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 " )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DNO_XPM -DNO_GENERAL -D_SEM_SEMUN_UNDEFINED")
add_executable(paradisgen ${PARADISGEN_SRCS})
add_executable(paradisrepart ${PARADISREPART_SRCS})
add_definitions(-DNO_XWINDOW )


target_link_libraries(paradisgen m geos pthread )
Install(TARGETS paradisgen
        RUNTIME DESTINATION ${PARADIS_SOURCE_DIR}/bin)

if(0)
message(STATUS "CMAKE_C_FLAGS = " ${CMAKE_C_FLAGS})
message(STATUS "CMAKE_C_FLAGS_DEBUG = " ${CMAKE_C_FLAGS_DEBUG})
message(STATUS "CMAKE_C_FLAGS_MINSIZEREL = " ${CMAKE_C_FLAGS_MINSIZEREL})
message(STATUS "CMAKE_C_FLAGS_RELEASE = " ${CMAKE_C_FLAGS_RELEASE})
message(STATUS "CMAKE_C_FLAGS_RELWITHDEBINFO = " ${CMAKE_C_FLAGS_RELWITHDEBINFO})

message(STATUS "CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
message(STATUS "CMAKE_CXX_FLAGS_DEBUG = " ${CMAKE_CXX_FLAGS_DEBUG})
message(STATUS "CMAKE_CXX_FLAGS_MINSIZEREL = " ${CMAKE_CXX_FLAGS_MINSIZEREL})
message(STATUS "CMAKE_CXX_FLAGS_RELEASE = " ${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO = " ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})

message(STATUS "CMAKE_EXE_LINKER_FLAGS = " ${CMAKE_EXE_LINKER_FLAGS})
message(STATUS "CMAKE_EXE_LINKER_FLAGS_DEBUG = " ${CMAKE_EXE_LINKER_FLAGS_DEBUG})
message(STATUS "CMAKE_EXE_LINKER_FLAGS_MINSIZEREL = " ${CMAKE_EXE_LINKER_FLAGS_MINSIZEREL})
message(STATUS "CMAKE_EXE_LINKER_FLAGS_RELEASE = " ${CMAKE_EXE_LINKER_FLAGS_RELEASE})
message(STATUS "CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO = " ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO})

message(STATUS "CMAKE_MODULE_LINKER_FLAGS = " ${CMAKE_MODULE_LINKER_FLAGS})
message(STATUS "CMAKE_MODULE_LINKER_FLAGS_DEBUG = " ${CMAKE_MODULE_LINKER_FLAGS_DEBUG})
message(STATUS "CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL = " ${CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL})
message(STATUS "CMAKE_MODULE_LINKER_FLAGS_RELEASE = " ${CMAKE_MODULE_LINKER_FLAGS_RELEASE})
message(STATUS "CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO = " ${CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO})

message(STATUS "CMAKE_SHARED_LINKER_FLAGS = " ${CMAKE_SHARED_LINKER_FLAGS})
message(STATUS "CMAKE_SHARED_LINKER_FLAGS_DEBUG = " ${CMAKE_SHARED_LINKER_FLAGS_DEBUG})
message(STATUS "CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL = " ${CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL})
message(STATUS "CMAKE_SHARED_LINKER_FLAGS_RELEASE = " ${CMAKE_SHARED_LINKER_FLAGS_RELEASE})
message(STATUS "CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO = " ${CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO})

message(STATUS "CMAKE_STATIC_LINKER_FLAGS = " ${CMAKE_STATIC_LINKER_FLAGS})
message(STATUS "CMAKE_STATIC_LINKER_FLAGS_DEBUG = " ${CMAKE_STATIC_LINKER_FLAGS_DEBUG})
message(STATUS "CMAKE_STATIC_LINKER_FLAGS_MINSIZEREL = " ${CMAKE_STATIC_LINKER_FLAGS_MINSIZEREL})
message(STATUS "CMAKE_STATIC_LINKER_FLAGS_RELEASE = " ${CMAKE_STATIC_LINKER_FLAGS_RELEASE})
message(STATUS "CMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO = " ${CMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO})

endif()
