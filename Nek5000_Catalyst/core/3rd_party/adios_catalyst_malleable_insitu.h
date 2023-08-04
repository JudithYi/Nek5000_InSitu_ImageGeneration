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


extern "C"{
    int adios_catalyst(MPI_Comm & newcomm, MPI_Comm & worldComm);
}

#endif
// HeaderTest-Exclude: adios_catalyst_melleable_inistu.h