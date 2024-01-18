#The Gamma Routing Flow Model

The Gamma Routing flow model is a concetpual model for flow propagation. 
It has been designed to work for complexe networks and can be easily coupled with distributed or semi-distributed hydrological model. 
Is run extremly fast. Thanks to the Gamma function which design the unit hydrogram, the model is dynamics and is able to reproduce complexe hydraulics behaviors.

The Gamma Routing model is differentiable and come with it adjoint. The code integrates a controleur and an optimizer (lbfgsb) to calibrate the parameters with a variationnal algorithm.

The Gamma Routing model is written in Fortran and come with a full Python interface. It can run and integrate Fortran or Python program.

Compilation:
make
n
y


