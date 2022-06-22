
#############################################################################
#
#    makefile: controls the build of various ParaDiS support utilities
#
#    Builds the following utilities:
#
#        paradisgen     --  problem generator
#        paradisrepart  --  creates a new problem decomposition with a new
#                           domain geometry from an previous nodal data file
#        paradisconvert --  Converts older format restart files to the
#                           current format
#        ctablegen      --  Serial version of utility used to create the
#                           FMM correction tables.
#        calcdensity    --  Calculates a dislocation density 'grid' from
#                           a specified restart file and writes it to a
#                           file for visualization via an external tool
#
#
#	NOTE: The utilities use various source modules from the parallel
#             code located in the parent directory.  These modules are
#             compiled in the parent directory with a different set of
#             preprocessor definitions than are needed here, so we create
#             links in this directory back to the source modlues and create
#             separate object modules for these sources.
#
############################################################################
include ../makefile.sys
include ../makefile.setup
include ../makefile.srcs

SRCDIR = ../src
INCDIR = ../include
BINDIR = ../bin


#
#	Define the exectutable, source and object modules for
#	the problem generator.
#

PARADISGEN     = paradisgen
PARADISGEN_BIN = $(BINDIR)/$(PARADISGEN)

PARADISGEN_CPP_SRCS = ParadisGen.cpp   \
               CreateConfig.cpp      \
               Decomp.cpp            \
               DLBfreeOld.cpp        \
               FindPreciseGlidePlane.cpp  \
               FMComm.cpp            \
               FMSigma2.cpp          \
               FMSupport.cpp         \
               FreeInitArrays.cpp    \
               GetNewNativeNode.cpp  \
               Heap.cpp              \
               InitCellDomains.cpp   \
               InitCellNatives.cpp   \
               InitCellNeighbors.cpp \
               InitSendDomains.cpp   \
               InitRemesh.cpp        \
               InitRemoteDomains.cpp \
               Matrix.cpp            \
               MemCheck.cpp          \
               Meminfo.cpp           \
               Param.cpp             \
               Parse.cpp             \
               PickScrewGlidePlane.cpp \
               QueueOps.cpp          \
               RBDecomp.cpp          \
               ReadRestart.cpp       \
               ReadBinaryRestart.cpp \
               RSDecomp.cpp          \
               Timer.cpp             \
               Util.cpp              \
      		   Math.cpp	           \
			   Grain.cpp

#PARADISGEN_CPP_SRCS = DisplayC.cpp   \
                   display.cpp 

PARADISGEN_SRCS = $(PARADISGEN_C_SRCS) $(PARADISGEN_CPP_SRCS)
PARADISGEN_OBJS = $(PARADISGEN_C_SRCS:.c=.o)  $(PARADISGEN_CPP_SRCS:.cpp=.o)


#
#	Define the exectutable, source and object modules for
#	the repartitioning utility.
#

PARADISREPART     = paradisrepart
PARADISREPART_BIN = $(BINDIR)/$(PARADISREPART)

PARADISREPART_CPP_SRCS = ParadisRepart.cpp  \
                    Decomp.cpp            \
                    DLBfreeOld.cpp        \
                    FindPreciseGlidePlane.cpp  \
                    FMComm.cpp            \
                    FMSigma2.cpp          \
                    FMSupport.cpp         \
                    FreeInitArrays.cpp    \
                    Heap.cpp              \
                    InitCellDomains.cpp   \
                    InitCellNatives.cpp   \
                    InitCellNeighbors.cpp \
                    InitRemoteDomains.cpp \
                    InitSendDomains.cpp   \
                    Matrix.cpp            \
                    MemCheck.cpp          \
                    Param.cpp             \
                    Parse.cpp             \
                    PickScrewGlidePlane.cpp \
                    QueueOps.cpp          \
                    ReadRestart.cpp       \
                    ReadBinaryRestart.cpp \
                    RBDecomp.cpp          \
                    RSDecomp.cpp          \
                    Timer.cpp             \
                    Util.cpp              \
                    WriteRestart.cpp		\
      		   		Math.cpp	           \
					Grain.cpp

#PARADISREPART_CPP_SRCS = DisplayC.cpp     \
                      display.cpp 

PARADISREPART_SRCS = $(PARADISREPART_C_SRCS) $(PARADISREPART_CPP_SRCS)
PARADISREPART_OBJS = $(PARADISREPART_C_SRCS:.c=.o)  $(PARADISREPART_CPP_SRCS:.cpp=.o)


#
#	Define the exectutable, source and object modules for
#	the density grid utility.
#

CALCDENSITY     = calcdensity
CALCDENSITY_BIN = $(BINDIR)/$(CALCDENSITY)

CALCDENSITY_CPP_SRCS = CalcDensity.cpp      \
                    Decomp.cpp            \
                    DLBfreeOld.cpp        \
                    FindPreciseGlidePlane.cpp  \
                    FMComm.cpp            \
                    FMSigma2.cpp          \
                    FMSupport.cpp         \
                    FreeInitArrays.cpp    \
                    Heap.cpp              \
                    InitCellDomains.cpp   \
                    InitCellNatives.cpp   \
                    InitCellNeighbors.cpp \
                    InitRemoteDomains.cpp \
                    InitSendDomains.cpp   \
                    Matrix.cpp            \
                    MemCheck.cpp          \
                    Param.cpp             \
                    Parse.cpp             \
                    PickScrewGlidePlane.cpp \
                    QueueOps.cpp          \
                    ReadRestart.cpp       \
                    ReadBinaryRestart.cpp \
                    RBDecomp.cpp          \
                    RSDecomp.cpp          \
                    Timer.cpp             \
                    Util.cpp              \
      		   		Math.cpp	           \
					Grain.cpp
             
#
#CALCDENSITY_CPP_SRCS = DisplayC.cpp     \
                       display.cpp 

CALCDENSITY_SRCS = $(CALCDENSITY_C_SRCS) $(CALCDENSITY_CPP_SRCS)
CALCDENSITY_OBJS = $(CALCDENSITY_C_SRCS:.c=.o)  $(CALCDENSITY_CPP_SRCS:.cpp=.o)


#
#	Define the exectutable for the serial version of the FMM
#       image correction table generator.  
#
#       NOTE: Except for the module containing main() for the ctablegen
#             utility, the source and object modules for this utility are
#             defined in ../makefile.srcs since they are needed in the
#             makefile in the main source directory to compile the parallel
#             version of the utility.
#

CTABLEGEN          = ctablegen
CTABLEGEN_BIN      = $(BINDIR)/$(CTABLEGEN)
CTABLEGEN_MAIN_SRC = CTableGen.cpp


#
#       Define the exectutable, source and object modules for
#       the generator for the stress tables used if fast-multipole
#       code is not used.
#

STRESSTABLEGEN     = stresstablegen
STRESSTABLEGEN_BIN = $(BINDIR)/$(STRESSTABLEGEN)

STRESSTABLEGEN_SRCS = StressTableGen.cpp \
                      FindPreciseGlidePlane.cpp \
                      Heap.cpp           \
                      Matrix.cpp         \
                      MemCheck.cpp       \
                      PickScrewGlidePlane.cpp \
                      QueueOps.cpp       \
                      Util.cpp      		\
      		   		  Math.cpp	           \
					  Grain.cpp

STRESSTABLEGEN_OBJS =  $(STRESSTABLEGEN_SRCS:.cpp=.o)


#
#       Define the exectutable, source and object modules for
#       the control parameter and data file converter (old style 
#       files to the current style)
#

PARADISCONVERT     = paradisconvert
PARADISCONVERT_BIN = $(BINDIR)/$(PARADISCONVERT)

PARADISCONVERT_CPP_SRCS = ParadisConvert.cpp \
                     Decomp.cpp            \
                     DLBfreeOld.cpp        \
                     FindPreciseGlidePlane.cpp  \
                     FMComm.cpp            \
                     FMSigma2.cpp          \
                     FMSupport.cpp         \
                     FreeInitArrays.cpp    \
                     Heap.cpp              \
                     InitCellDomains.cpp   \
                     InitCellNatives.cpp   \
                     InitCellNeighbors.cpp \
                     InitRemoteDomains.cpp \
					 InitSendDomains.cpp   \
                     Matrix.cpp            \
                     MemCheck.cpp          \
                     Param.cpp             \
                     Parse.cpp             \
                     PickScrewGlidePlane.cpp \
                     QueueOps.cpp          \
                     ReadRestart.cpp       \
                     ReadBinaryRestart.cpp \
                     RBDecomp.cpp          \
                     RSDecomp.cpp          \
                     Timer.cpp             \
                     Util.cpp         	\
                     Math.cpp         	\
					 Grain.cpp

#PARADISCONVERT_CPP_SRCS = DisplayC.cpp     \
                       display.cpp 

PARADISCONVERT_SRCS = $(PARADISCONVERT_C_SRCS) $(PARADISCONVERT_CPP_SRCS)
PARADISCONVERT_OBJS = $(PARADISCONVERT_C_SRCS:.c=.o) $(PARADISCONVERT_CPP_SRCS:.cpp=.o)


###########################################################################
#
#	Define a rule for converting .c files to object modules.
#	All modules are compile serially in this directory
#
###########################################################################

.c.o:		makefile ../makefile.sys ../makefile.setup
		$(CC_SERIAL) $(OPT) $(CCFLAG_SERIAL) $(INCS_SERIAL) -c $<

.C.o:		makefile ../makefile.sys ../makefile.setup
		$(CPP_SERIAL) $(OPT) $(CPPFLAG_SERIAL) $(INCS_SERIAL) -c $<


###########################################################################
#
#	Define all targets and dependencies below
#
###########################################################################

all:		$(PARADIS_SRCS) $(CTABLEGEN_SRCS) $(CTABLEGEN_MAIN_SRC) $(CTABLEGEN) $(PARADISGEN) \
			$(PARADISREPART) $(STRESSTABLEGEN) $(PARADISCONVERT) \
			$(CALCDENSITY)

clean:
		rm -f *.o $(PARADIS_SRCS) $(CTABLEGEN_SRCS) $(CTABLEGEN_MAIN_SRC) $(CTABLEGEN_BIN) \
			$(PARADISGEN_BIN) $(PARADISCONVERT_BIN) \
			$(PARADISREPART_BIN) $(STRESSTABLEGEN_BIN) \
			$(CALCDENSITY_BIN)

depend:		 *.c $(SRCDIR)/*.c $(INCDIR)/*.h makefile
		makedepend -Y$(INCDIR) *.c  -fmakefile.dep

$(BINDIR):
		mkdir $(BINDIR)

#
#	Create any necessary links in the current directory to source
#	modules located in the SRCDIR directory
#

$(PARADIS_SRCS): $(SRCDIR)/$@
		- @ ln -s  -f $(SRCDIR)/$@ ./$@ > /dev/null 2>&1

$(CTABLEGEN_MAIN_SRC): $(SRCDIR)/$@
		- @ ln -s  -f $(SRCDIR)/$@ ./$@ > /dev/null 2>&1


#
#       Targets for each of the utility executables
#

$(PARADISCONVERT):	$(BINDIR) $(PARADISCONVERT_BIN)
$(PARADISCONVERT_BIN):	$(PARADISCONVERT_SRCS) $(PARADISCONVERT_OBJS)
		$(CPP_SERIAL) $(CPPFLAG_SERIAL) $(INCS_SERIAL) $(OPT) -o $@ \
			$(PARADISCONVERT_OBJS) $(LIB_SERIAL)

$(STRESSTABLEGEN):	$(BINDIR) $(STRESSTABLEGEN_BIN)
$(STRESSTABLEGEN_BIN):	$(STRESSTABLEGEN_SRCS) $(STRESSTABLEGEN_OBJS)
		$(CPP_SERIAL) $(CCFLAG_SERIAL) $(INCS_SERIAL) $(OPT) -o $@ \
		$(STRESSTABLEGEN_OBJS) $(LIB_SERIAL)

$(CTABLEGEN):	$(BINDIR) $(CTABLEGEN_BIN)
$(CTABLEGEN_BIN):	$(CTABLEGEN_SRCS) $(CTABLEGEN_OBJS)
		$(CPP_SERIAL) $(CPPFLAG_SERIAL) $(INCS_SERIAL) $(OPT) -o $@ \
		$(CTABLEGEN_OBJS) $(LIB_SERIAL)

$(PARADISGEN):	$(BINDIR) $(PARADISGEN_BIN)
$(PARADISGEN_BIN):	$(PARADISGEN_SRCS) $(PARADISGEN_OBJS)
		$(CPP_SERIAL) $(CPPFLAG_SERIAL) $(INCS_SERIAL) $(OPT) -o $@ \
		$(PARADISGEN_OBJS) $(LIB_SERIAL)

$(PARADISREPART):	$(BINDIR) $(PARADISREPART_BIN)
$(PARADISREPART_BIN):	$(PARADISREPART_SRCS) $(PARADISREPART_OBJS)
		$(CPP_SERIAL) $(CPPFLAG_SERIAL) $(INCS_SERIAL) $(OPT) -o $@ \
		$(PARADISREPART_OBJS) $(LIB_SERIAL)

$(CALCDENSITY):	$(BINDIR) $(CALCDENSITY_BIN)
$(CALCDENSITY_BIN):	$(CALCDENSITY_SRCS) $(CALCDENSITY_OBJS)
		$(CPP_SERIAL) $(CPPFLAG_SERIAL) $(INCS_SERIAL) $(OPT) -o $@ \
		$(CALCDENSITY_OBJS) $(LIB_SERIAL)
