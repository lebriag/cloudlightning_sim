# CloudLightning Simulator

CloudLightning Simulation is a generalized and extensible simulation framework that enables the seamless simulation and experimentation of emerging Cloud computing infrastructures and  HPC applications. The framework is inherently parallel; it is written using the C/C++ programming language and the MPI and OpenMP APIs and enables the exploitation of distributed and shared memory parallel techniques for the acceleration of Cloud simulation.

Instead of following a Discrete Events approach as other Cloud simulation frameworks, the proposed work is based on a time-advancing loop where the status of the Cloud system can change at each time step. Using this approach allows for reduced memory requirements, since events need not be created or stored; instead the system reacts to incoming tasks based on a prescribed time interval.

The project documentation can be found online at:

     http://www.iti.gr/~kouzinopoulos/

## Getting started

Clone the project:

     https://bitbucket.org/cloudlightning/cloudlightning-simulator

## Prerequisites

* Required gcc versions 4.8-6.0 (https://gcc.gnu.org/)
* Required OpenMPI version 2.1 (https://www.open-mpi.org)

Install OpenMPI

     $ ./configure --prefix=/where/to/install
     $ make all install
     $ sudo apt-get install libopenmpi-dev

## Input and Output Data

Input JSON data files are located at the corresponding *input* folder of the project. In *input* folder there are three JSON files containing the appropriate parameters about
applications, brokers and cells. Each of these parameters can be configured by the user before running the simulation. *CellData.json* file include some global information such
as the maximum time of the simulation, the update interval of the simulation, the number of the cells and the integration (or not) of the SOSM system.

In *output* folder will be generated the *outputCLsim.json* file after running the simulator.

##Project Build

Navigate to the directory of the project, open the terminal and enter

     $ cmake .
     $ make

##Documentation Build (optional)

To build the project documentation locally, navigate to the directory of the project, open the terminal and enter

     $ make doc

Then, open the doc/html/index.html file.

##Run the CL Simulator

At the directory of the project, open the terminal and enter

     $ sh cl_sim.out

After running the simulation, an *outputCLsim.json* file will be created in the *output* folder. This file can be used as an input for the GUI implementation of the CL Simulator in order
to visualize the results of the simulator.

The Simulator Visualization tool can be found at the following url:

     https://bitbucket.org/cloudlightning/cl-simulatorvisualization

##Deliverable

The deliverable 7.1.1 of the CloudLightning Project can be found at the following url:

     http://cloudlightning.eu/work-packages/public-deliverables/
