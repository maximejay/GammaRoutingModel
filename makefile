#GammaRouting is a conceptual flow propagation model
#Copyright 2022 Maxime Jay-Allemand

#This file is part of GammaRouting.

#GammaRouting is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

#GammaRouting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with GammaRouting. If not, see <https://www.gnu.org/licenses/>

#Contact: maxime.jay.allemand@hydris-hydrologie.fr



#Compiler and Linker
CC := gcc

#The Target Binary Program
TARGET     := routing

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIRc     := src/c
SRCDIRf90   := src/f90
SRCDIRf77   := src/f77
SRCDIRwrap  := src/py/wrappers
SRCDIRwrapping  := src/py/gamma/wrapping
SRCDIRpy    := src/py
INCDIR      := .
BUILDDIR    := obj
BUILDDIRC   := obj
BUILDDIRF90 := obj
BUILDDIRF77 := obj
TARGETDIR   := bin
RESDIR      := res
SRCEXTC     := c
SRCEXTF90   := f90
SRCEXTF77   := f
DEPEXT      := d
OBJEXT      := o



#Flags, Libraries and Includes
#CFLAGS      := -g -O3 -march=native -ffast-math


#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
SOURCESC    := $(shell find $(SRCDIRc) -type f -name "*.$(SRCEXTC)")
OBJECTSC    := $(patsubst $(SRCDIRc)/%,$(BUILDDIRC)/%,$(SOURCESc:.$(SRCEXTC)=.$(OBJEXT)))

SOURCESF90  := $(shell find $(SRCDIRf90) -type f -name "*.$(SRCEXTF90)")
OBJECTSF90  := $(patsubst $(SRCDIRf90)/%,$(BUILDDIRF90)/%,$(SOURCESF90:.$(SRCEXTF90)=.$(OBJEXT)))

SOURCESF77  := $(shell find $(SRCDIRf77) -type f -name "*.$(SRCEXTF77)")
OBJECTSF77  := $(patsubst $(SRCDIRf77)/%,$(BUILDDIRF77)/%,$(SOURCESF77:.$(SRCEXTF77)=.$(OBJEXT)))

F90_INTERFACE := src/f90/mod_routing_setup.f90 src/f90/mod_routing_mesh.f90 src/f90/mod_routing_states.f90 src/f90/mod_routing_results.f90 src/f90/mod_routing_parameters.f90 src/f90/mod_gamma_interface.f90

OBJ_INTERFACE := $(BUILDDIR)/*.o 
F90_WRAPPERS  := $(SRCDIRwrap)/*.f90
SHAREDLIB     := wrapping
KIND_MAP      := kind_map
DOC_PLUGIN    := doc_plugin.py


#Check the compiler version and add extra flag if needed
GCCVERSION = $(shell gcc -dumpversion)
GCC_VERSION_GREATER_THAN_9 = $(shell expr `gcc -dumpversion | cut -d . -f 1` \> 9)
ADDF90FLAGS = 
ifeq "$(GCC_VERSION_GREATER_THAN_9)" "1"
	ADDF90FLAGS += -fallow-argument-mismatch
	ADDF77FLAGS += -fallow-argument-mismatch
endif

#Default Make
all:
	@echo
	@echo "BUILD OPTIONS"
	@echo "============="
	@echo
	@echo "Using Gfortran version "$(GCCVERSION)
	@echo "Extra compiler flags = "$(ADDF90FLAGS)
	@echo
	
	@echo "Would you like to regenerate the adjoint? [Y/n]" ; \
	echo "   Default: n" ; \
	read ans_adj ; \
	echo "" ; \
	echo "Would you like to generate the Python interface? [Y/n]" ; \
	echo "   Default: n" ; \
	read ans_python ; \
	if [ $${ans_adj} = "y" ]; then \
	  make adjoint ; \
	elif [ $${ans_adj} = "Y" ]; then \
	  make adjoint ; fi ; \
	if [ $${ans_python} = "y" ]; then \
	  make f90py ; \
	elif [ $${ans_python} = "Y" ]; then \
	  make f90py ; \
	else make f90 ; fi ;


#Compile Fortran
f90: directories target_gnu

#Compile Fortran and Python
f90py: f90 wrappers module package

#Compile Fortran debug mode
f90debug: directories target_debug

#Compile Fortran debug mode and Python
f90pydebug: f90debug wrappers module package


target_gnu: FC = gfortran
# http://www.fortran90.org/src/faq.html#what-compiler-options-should-i-use-for-development
target_gnu: F90FLAGS      = -cpp -O3 -march=native -funroll-loops -fPIC $(ADDF90FLAGS)
target_gnu: F77FLAGS      = -O3 -march=native -funroll-loops -fPIC $(ADDF90FLAGS)
target_gnu: CFLAGS        := -g -O3 -march=native -fPIC
target_gnu: directories $(TARGET) 


wrappers:
	@echo "********************************************"
	@echo ""
	@echo " Making wrappers "
	@echo ""
	@echo "********************************************"
	f90wrap -m $(SHAREDLIB) $(F90_INTERFACE) -k $(KIND_MAP) --documentation-plugin $(DOC_PLUGIN) --package
	mv f90wrap_*.f90 $(SRCDIRwrap)/.

module:
	@echo "********************************************"
	@echo ""
	@echo " Making module extension "
	@echo ""
	@echo "********************************************"
	f2py-f90wrap -c --fcompiler=gfortran --f90flags='-cpp -fPIC -fmax-errors=1 -Iobj' --arch='-march=native' --opt='-O3 -funroll-loops' --build-dir . -m _$(SHAREDLIB) $(F90_WRAPPERS) $(F90_INTERFACE)
	mv $(SHAREDLIB)/mod* $(SRCDIRwrapping)/.
	mv _$(SHAREDLIB)* $(SRCDIRwrapping)/.
	rm -rf $(SHAREDLIB)

package:
	@echo "********************************************"
	@echo ""
	@echo " Installing Python smash package "
	@echo ""
	@echo "********************************************"
	pip install --no-compile src/py/.
	
package_edit:
	@echo "********************************************"
	@echo ""
	@echo " Installing Python smash package (Editable) "
	@echo ""
	@echo "********************************************"
	pip install -e src/py/.


#Mode debug - set appropriate compiler flags : -Wimplicit-interface
target_debug: FC = gfortran
target_debug: F90FLAGS      = -Wall -Wextra -fPIC -fmax-errors=1 -cpp -g -fcheck=all -fbacktrace $(ADDF90FLAGS)
target_debug: F77FLAGS      = -Wall -Wextra -fPIC -fmax-errors=1 -cpp -g -fcheck=all -fbacktrace $(ADDF90FLAGS)
target_debug: CFLAGS        = -g -fbacktrace -fPIC
target_debug: directories $(TARGET) 

tests:
	cd tests
	
#Remake
remake: cleaner all

#Make the Directories
directories:
	@echo "********************************************"
	@echo ""
	@echo " Making directories "
	@echo ""
	@echo "********************************************"
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(SRCDIRwrap)
	@mkdir -p $(SRCDIRwrapping)


#Clean only Objecst
clean:
	@$(RM) -rf $(BUILDDIR)
	@$(RM) -rf *.mod

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) -rf $(TARGETDIR)
	@$(RM) -rf *.mod
	@$(RM) -rf $(SRCDIRwrapping)/mod*
	@$(RM) -rf $(SRCDIRwrapping)/_$(SHAREDLIB)*
	@$(RM) -rf $(SRCDIRwrap)


adjoint:
	@echo "********************************************"
	@echo ""
	@echo " Making adjoint "
	@echo ""
	@echo "********************************************"
	cd differentiation/; ./make_adjoint.sh


#Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
$(TARGET): \
 obj/adStack.o \
 obj/adBuffer.o \
 obj/mod_gamma_sorting.o \
 obj/mod_routing_setup.o \
 obj/mod_routing_mesh.o \
 obj/mod_routing_parameters.o \
 obj/mod_routing_states.o \
 obj/mod_routing_results.o \
 obj/mod_gamma_function.o \
 obj/mod_gamma_routing.o \
 obj/mod_gamma_interface.o \
 obj/cost_function.o \
 obj/run_forward.o \
 obj/lbfgsb.o \
 obj/AADJ_b.o \
 obj/ATLM_d.o \
 obj/control.o \
 $(OBJECTSFC) \
 $(OBJECTSF77) \
 $(OBJECTSF90)
	$(FC) $(LDFLAGS) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

$(BUILDDIRF90)/%.$(OBJEXT): $(SRCDIRf90)/%.$(SRCEXTF90)
	@mkdir -p $(dir $@)
	$(FC) $(F90FLAGS) $(INC) -c -o $@ $<

$(BUILDDIRF77)/%.$(OBJEXT): $(SRCDIRf77)/%.$(SRCEXTF77)
	$(FC) $(F77FLAGS) -c -o $@ $<

$(BUILDDIRC)/%.$(OBJEXT): $(SRCDIRc)/%.$(SRCEXTC)
	$(CC) $(CFLAGS) -c -o $@ $<


#Non-File Targets
.PHONY: all remake clean cleaner resources fortrangis

