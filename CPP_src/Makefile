#***************************************************************************/
#**    Author: Joey Azofeifa                                               */
#**    Date:   1/20/2016	                                           */
#**    Description: Makefile for SE Calculator                             */
#***************************************************************************/

CXX             = mpic++
CXXFLAGS        = -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11 -fopenmp -Wno-write-strings -Wno-literal-suffix 
EXEC            = ${PWD}/SE
ARCH  		= getconf LONG_BIT

CPP_FLAGS_32  	= -D32_BIT
CPP_FLAGS_64  	= -D64_BIT

GCCVERSION 		= $(shell ${CXX} -dumpversion)
SE: main.o  simulate.o load.o split.o read_in_parameters.o get_motif_pvalues.o scanner.o out.o ACGT_profile.o MPI_comm.o collect_sample_statistics.o error_stdo_logging.o
	@printf "linking           : "
	@${CXX} ${CXXFLAGS} ${PWD}/collect_sample_statistics.o ${PWD}/simulate.o ${PWD}/scanner.o ${PWD}/ACGT_profile.o ${PWD}/MPI_comm.o  ${PWD}/main.o ${PWD}/load.o ${PWD}/out.o ${PWD}/split.o ${PWD}/read_in_parameters.o ${PWD}/error_stdo_logging.o ${PWD}/get_motif_pvalues.o -o ${EXEC}
	@printf "done\n"
	@echo "========================================="
	@printf "finished :)\n"




main.o: 
	@echo "========================================="
	@echo "GCC version: "${GCCVERSION}
	@printf "main              : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/main.cpp 
	@printf "done\n"
error_stdo_logging.o:
	@printf "stdo/logging      : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/error_stdo_logging.cpp
	@printf "done\n"

collect_sample_statistics.o:
	@printf "CSS               : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/collect_sample_statistics.cpp
	@printf "done\n"

MPI_comm.o:
	@printf "MPI_comm          : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/MPI_comm.cpp
	@printf "done\n"


simulate.o:
	@printf "simulate          : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/simulate.cpp 
	@printf "done\n"



ACGT_profile.o:
	@printf "ACGT_profile      : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/ACGT_profile.cpp 
	@printf "done\n"
scanner.o:
	@printf "scanner           : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/scanner.cpp 
	@printf "done\n"

out.o:
	@printf "out               : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/out.cpp 
	@printf "done\n"


get_motif_pvalues.o: 
	@printf "get_motif_pvalues : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/get_motif_pvalues.cpp 
	@printf "done\n"


load.o: 
	@printf "load              : "	
	@${CXX} -c ${CXXFLAGS} ${PWD}/load.cpp 
	@printf "done\n"
split.o:
	@printf "split             : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/split.cpp 
	@printf "done\n"


read_in_parameters.o:
	@printf "read_in_parameters: "
	@${CXX} -c ${CXXFLAGS} ${PWD}/read_in_parameters.cpp 
	@printf "done\n"





clean:
	@rm -f *.o