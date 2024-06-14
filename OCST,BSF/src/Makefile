################################################################################
# When you adapt this makefile to compile your CPLEX programs
# set CPLEXDIR to the directories where CPLEX is installed:
LIB = ../../../../../encs/pkg/cplex-22.1.0/root/cplex
OCST: main.o read_data.o output_data.o Heuristics.o Benders.o Benders_master.o Benders_subproblem.o CLM_utils.o CLM_graph.o CLM_OCSTP.o backwardAlgorithm.o FortificationProblem.o newSampling.o perceiveDamage.o PreFortification.o realDamage.o sampling.o backwardAlgorithmfinal.o storeWaitList.o addCoverconstraint.o
	gcc -O2 -I$(LIB)/include main.o read_data.o output_data.o Heuristics.o Benders.o Benders_master.o Benders_subproblem.o CLM_utils.o CLM_graph.o CLM_OCSTP.o backwardAlgorithm.o FortificationProblem.o newSampling.o perceiveDamage.o PreFortification.o realDamage.o sampling.o backwardAlgorithmfinal.o storeWaitList.o addCoverconstraint.o -L$(LIB)/lib/x86-64_linux/static_pic  -lcplex -lm -lpthread -ldl -o OCST

main.o: main.c headers.h
	gcc -O2 -c main.c
read_data.o: read_data.c headers.h
	gcc -O2 -c read_data.c
output_data.o: output_data.c headers.h
	gcc -O2 -c output_data.c
Heuristics.o: Heuristics.c headers.h
	gcc -O2 -c Heuristics.c
Benders.o: Benders.c headers.h
	gcc -O2 -c Benders.c
Benders_master.o: Benders_master.c headers.h
	gcc -O2 -c Benders_master.c
Benders_subproblem.o: Benders_subproblem.c headers.h
	gcc -O2 -c Benders_subproblem.c 
CLM_utils.o: CLM_utils.c headers.h
	gcc -O2 -c CLM_utils.c
CLM_graph.o: CLM_graph.c headers.h
	gcc -O2 -c CLM_graph.c
CLM_OCSTP.o: CLM_OCSTP.c headers.h
	gcc -O2 -c CLM_OCSTP.c
backwardAlgorithm.o: backwardAlgorithm.c headers.h
	gcc -O2 -c backwardAlgorithm.c
FortificationProblem.o: FortificationProblem.c headers.h
	gcc -O2 -c FortificationProblem.c
newSampling.o: newSampling.c headers.h
	gcc -O2 -c newSampling.c
perceiveDamage.o: perceiveDamage.c headers.h
	gcc -O2 -c perceiveDamage.c
PreFortification.o: PreFortification.c headers.h
	gcc -O2 -c PreFortification.c
realDamage.o: realDamage.c headers.h
	gcc -O2 -c realDamage.c
sampling.o: sampling.c headers.h
	gcc -O2 -c sampling.c
backwardAlgorithmfinal.o: backwardAlgorithmfinal.c headers.h
	gcc -O2 -c backwardAlgorithmfinal.c
addCoverconstraint.o: addCoverconstraint.c headers.h
	gcc -O2 -c addCoverconstraint.c
storeWaitList.o: storeWaitList.c headers.h
	gcc -O2 -c storeWaitList.c
clean :
	rm MST main.o read_data.o output_data.o Heuristics.o Benders.o Benders_master.o Benders_subproblem.o CLM_utils.o CLM_graph.o CLM_OCSTP.o backwardAlgorithm.o FortificationProblem.o newSampling.o perceiveDamage.o PreFortification.o realDamage.o sampling.o backwardAlgorithmfinal.o storeWaitList.o addCoverconstraint.o
