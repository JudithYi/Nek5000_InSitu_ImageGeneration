# How to try Nek5000 with in-situ image generation? 
1. Install Paraview 5.9 with Catalyst with the instructions (3.1-3.3) on: 
https://github.com/KTH-Nek5000/InSituPackage

2. Compile and run the synchronous in-situ image generation with the code in folder "sync":
2.0 Copy the 'tgv' folder in 'Nek5000_Catalyst/examples/' folder into 'Nek5000_Catalyst/run' folder and copy all the files in 'sync' folder into 'Nek5000_Catalyst/run/tgv' folder where codes need to be compiled and run.
2.1 Change the path in make.sh and run.sh
2.2 Compile with make.sh 
2.3 Copy the python scripts in pyfile in 'pyfile' folder into 'Nek5000_Catalyst/run/tgv' folder.
2.4 Run the code with run.sh after substituting tgv.par in 'tgv' folder with the tgv.par in this folder.
2.5 Try the different number of cores for nek5000 (in run.sh) 

3. Compile and run the asynchronous in-situ image generation with the code in folder "async"
3.0 Copy the 'tgv' folder in 'Nek5000_Catalyst/examples/' folder into 'Nek5000_Catalyst/run' folder and copy all the files in 'async' folder into 'Nek5000_Catalyst/run/tgv' folder where codes need to be compiled and run.
3.1 Change the path in make.sh (also the one in "insitu" folder) and run.sh
3.2 Compile with make.sh (maybe also make.sh in "insitu" folder) and copy the catalystSpace into 'Nek5000_Catalyst/run/tgv' folder.
3.3 Copy the python scripts in pyfile in 'pyfile' folder into 'Nek5000_Catalyst/run/tgv' folder.
3.4 Run the code with run.sh after substituting tgv.par in 'Nek5000_Catalyst/run/tgv' folder with the tgv.par in this folder. 
3.5 Try the different number of cores for nek5000 and catalystSpace (in run.sh), engine type, other parameters (in config/config.xml) with the help of the information from ADIOS website


File changed compared to the original Nek5000:
1. core/drive.f <- the main function with split MPI communicator
2. core/drive1.f <- deactivated prepost function
3. core/3rd_party/nek_in_situ.f <- additional in-situ functions
4. core/3rd_party/nek_catalyst.cxx + catalyst.f <- additional files for synchronous in-situ image generation. 
5. core/3rd_party/nek_adios.cpp + adios.f <- additional files for asynchronous in-situ image generation. 
6. core/makefile.template + makenek.inc <- modified compile files 