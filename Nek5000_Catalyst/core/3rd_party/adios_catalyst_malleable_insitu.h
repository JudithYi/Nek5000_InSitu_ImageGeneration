#ifndef adios_catalyst_melleable_inistu_h
#define adios_catalyst_melleable_inistu_h

#include <adios2.h>
#include <mpi.h>
#include <string>
#include <vector>
#include <iostream>
#include <time.h>
#include <algorithm>
#include "nek_catalyst_async.h"
#include "myCPPythonAdaptorAPI.h"

int adios_catalyst_init(
    MPI_Comm & comm_in, 
    MPI_Comm & worldComm, 
    std::string enginePair, 
    const int firstPair, 
    const int group=0
);
int adios_catalyst();
//int adios_catalyst_run();

#endif
// HeaderTest-Exclude: adios_catalyst_melleable_inistu.h
