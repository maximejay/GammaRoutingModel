#!/bin/bash
#Parse input arguments

DEBUG=0
F90BUILD=0
PYWRAPPERS=0
SHAREDLIB=0
INSTALL=0
CLEAN=0
EDIT=0

i=1
until [ $# = 0 ]
do
case $1 in
    --all)
    shift
    F90BUILD=1
    PYWRAPPERS=1
    SHAREDLIB=1
    INSTALL=1
    CLEAN=1
    ;;
    --debug)
    shift
    DEBUG=1
    ;;
    --edit)
    shift
    EDIT=1
    ;;
    --f90build)
    shift
    F90BUILD=1
    ;;
    --pywrappers)
    shift
    PYWRAPPERS=1
    ;;
    --sharedlib)
    shift
    SHAREDLIB=1
    ;;
    --install)
    shift
    INSTALL=1
    ;;
    --clean)
    shift
    CLEAN=1
    ;;
    *)
    echo "Unknow option $1 ..."
    echo "Use --help to see availbale options."
    exit 0
    ;;
esac
done



#make the adjoint
#~ cd differentiation/; ./make_adjoint.sh

#directories
src_dir=$(pwd)/src/
src_f90=$src_dir/f90/
src_f77=$src_dir/f77/
src_c=$src_dir/c/
src_py=$src_dir/py/
src_f90wrap=$src_dir/f90wrap/

wrapped_lib_name="libfgamma"
build_dir=$(pwd)/build/

if [ $CLEAN -eq 1 ] ; then
    rm -r $src_py/gamma_routing_model/$wrapped_lib_name
    rm -r $build_dir
fi

mkdir -p $build_dir
cd $build_dir
#creation d'un module python
mkdir -p $wrapped_lib_name

cp $src_f90/*.f90 .
cp $src_f77*.f .
cp $src_c/*.c .
cp $src_f90wrap/kind_map .

#all source files being compiled: order matters
cfiles="adStack.c"
f77files="lbfgsb.f adBuffer.f"
f90files="mod_routing_setup.f90 mod_routing_mesh.f90 mod_routing_parameters.f90 mod_routing_states.f90 mod_routing_results.f90 mod_gamma_function.f90 mod_gamma_routing.f90 mod_gamma_sorting.f90 mod_gamma_interface.f90 generic_subroutine_gamma_routing.f90 generic_subroutine_gamma_function.f90 run_forward.f90 AADJ_b.f90 ATLM_d.f90 control.f90 cost_function.f90 gradient_test.f90"
f90mainprog="routing.f90"

#module fortran file that going to be wrapped with python: order matters
f90files_wrap="mod_routing_setup.f90 mod_routing_mesh.f90 mod_routing_parameters.f90 mod_routing_states.f90 mod_routing_results.f90 mod_gamma_interface.f90"


#compiler flags
if [ $DEBUG -eq 1 ] ; then
    c_flags="-fPIC -g -fbacktrace -fPIC"
    f_flags="-fPIC -Wall -Wextra -fPIC -fmax-errors=1 -cpp -g -fcheck=all -fbacktrace -fallow-argument-mismatch"
else
    c_flags="-fPIC -O3"
    f_flags="-fPIC -O3 -march=native -funroll-loops -fallow-argument-mismatch"
fi

if [ $PYWRAPPERS -eq 1 ] ; then
    #Make the python wrappers and create a python module libfgamma
    f90wrap -k "kind_map" -m "$wrapped_lib_name" $f90files_wrap
    
    #change how the lib is imported
    sed -i "s/^\(import \_${wrapped_lib_name}\)/from \. \1/g" ${wrapped_lib_name}.py
    mv ./${wrapped_lib_name}.py ${wrapped_lib_name}/.
fi

if [ $F90BUILD -eq 1 ] ; then
    #compilation des sources et créations des objets .o
    gcc $c_flags -c -g -fPIC $cfiles
    gfortran $f_flags -c $f77files
    gfortran $f_flags -c $f90files

    #compilation of the main prog
    gfortran $f90mainprog *.o -o $(basename f90mainprog)
fi

if [ $SHAREDLIB -eq 1 ] ; then

    #order matters 
    f90wrapped_files=""
    for file in $f90files_wrap
    do
        f90wrapped_files+=" f90wrap_${file}"
    done
    
    #Création de la librairie partagée
    FFLAGS="-I$(pwd)" FC=gfortran f2py -c --f90flags="$f_flags" -m "_${wrapped_lib_name}" $f90wrapped_files *.o
    
    #create the __init__.py module file
    echo "from .libfgamma import Mod_Gamma_Interface, Mod_Gamma_Routing_Setup, Mod_Gamma_Routing_Mesh, Mod_Gamma_Routing_Parameters, Mod_Gamma_Routing_Results, Mod_Gamma_Routing_States" > ${wrapped_lib_name}/__init__.py

    #move all component inside the module directories and move it inside the source module
    mv _${wrapped_lib_name}* ${wrapped_lib_name}/.
fi

#Before install the gamma module, copy the wrappers and the shared lib module in the source code
cp -u -r ./${wrapped_lib_name} $src_py/gamma_routing_model/.

if [ $INSTALL -eq 1 ] ; then
    cd $src_dir
    if [ $EDIT -eq 1 ] ; then
        #intall the module with editable mode
        pip install -e ./py
    else
        pip install ./py
    fi
fi
