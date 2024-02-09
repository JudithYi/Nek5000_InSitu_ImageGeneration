#ifndef nek_adios2_h
#define nek_adios2_h

#include <adios2.h>
#include <mpi.h>
#include <string>
#include <iostream>
#include <ctime>
#include <vector>

void init_multiple_type(const int type_num_in);

void adios_writer_init(
    MPI_Comm & comm,
    MPI_Comm & worldComm, 
    std::string engineName,
    const int iostep_in,
    const int type_index
);

#endif
