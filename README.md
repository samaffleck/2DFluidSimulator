# 2DFluidSimulator
2DFluidSimulator is a static libary for solving fluid flow PDEs using the finite volume method.

# Build 
Build system uses cmake
1. Clone the repositry: git clone 
2. Create the "build" directory, and change directory into "build"
3. Run: cmake ..

# Linking to the libary
Once you have built the libary, the .lib file will be in the build directory. You will need to include the "include" folder and link to the correct libary in your project

# Testing
In "tests/" there are test files which can be used as examples for how to set up a simulation.
