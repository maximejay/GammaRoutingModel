#Differentiation
#You need Tapenade to run this code : https://www-sop.inria.fr/tropics/tapenade/
#invoke this script by the commande : make all

#This script differentiate the following subroutines :
FILES_FTR= mod_routing_mesh.f90 mod_routing_setup.f90 mod_routing_parameters.f90 mod_routing_states.f90 mod_routing_results.f90 mod_gamma_routing.f90 cost_function.f90 run_forward.f90
#mod_gamma_function.f90
#Option -fixinterface solves the following warning : PROBLEMATIC INTERACTION WITH INDEPENDENT AND DEPENDENT PARAMETER SETS: https://www-sop.inria.fr/tropics/tapenade/faq.html#validation

#Path to a local copie of Tapenade
TAPENADE=./tapenade_3.16/

all: d b

d: $(FILES_FTR)
	tapenade -forward -fixinterface -msglevel 20 -o TLM -head routing_hydrogram_forward \
	-vars "routing_parameter" \
	-outvars "cost" \
	-O $(PWD)/forward $(FILES_FTR)
	
	
b: $(FILES_FTR)
	tapenade -backward -fixinterface -msglevel 20 -o ADJ -head routing_hydrogram_forward \
	-vars "routing_parameter" \
	-outvars "cost" \
	-O $(PWD)/backward $(FILES_FTR)
	
	
dJ_dQin_adj: $(FILES_FTR)
	tapenade -backward -fixinterface -msglevel 20 \
	-head "routing_hydrogram_forward(routing_parameter)\(cost) routing_hydrogram_forward(inflows)\(cost)" \
	-output "ADJ" \
	-O $(PWD)/backward $(FILES_FTR)

dJ_dQin_tlm: $(FILES_FTR)
	tapenade -forward -fixinterface -msglevel 20 \
	-head "routing_hydrogram_forward(routing_parameter)\(cost) routing_hydrogram_forward(inflows)\(cost)" \
	-output "TLM" \
	-O $(PWD)/forward $(FILES_FTR)
	
	
c: $(FILES_FTR)
	tapenade -p -o nodifferentiated $(FILES_FTR)
	

v: $(FILES_FTR)
	tapenade -forward -context -fixinterface -msglevel 20 -o ttest -head routing_hydrogram_forward \
	-vars "routing_parameter" \
	-outvars "cost" \
	-O $(PWD)/forward $(FILES_FTR)
	
	
clean:
	rm -f forward/*
	rm -f backward/*
	
