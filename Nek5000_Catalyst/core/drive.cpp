#include "mpi.h"
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <chrono>
#include <thread>
#include "3rd_party/adios_catalyst_malleable_insitu.h"
#include "3rd_party/nek_adios2.h"

/* Nek is wrtier here and programmed in fortran*/
extern "C"{
	void nek_init_malleable_insitu_(int * comm);
	void nek_solve_malleable_insitu_first_(int * numInsitu);
	void in_situ_init_(int * worldComm);
	void nek_solve_malleable_insitu_(int * numInsitu);
	void nek_end_();
}

int main(int argc, char* argv[]){

    MPI_Comm comm = MPI_COMM_NULL;
    MPI_Comm comm_reader = MPI_COMM_NULL;
    std::vector<MPI_Comm> comm_writer;
    
    std::string enginePair;
    int colour;
    std::vector<int> colours;
    int comm_f;
    int world_rank, world_size;
    int cores_per_socket = std::stoi(argv[1]);
    int cores_in_situ = std::stoi(argv[2]);
    int num_groups = std::stoi(argv[3]);
    colours.resize(num_groups);
    int i;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if(world_rank%cores_per_socket>=num_groups*cores_in_situ){
        colour = num_groups;
        for(i = 0; i < num_groups; ++i){
            colours[i] = 1;
        }
    }else{
        colour = static_cast<int>(world_rank%cores_per_socket/cores_in_situ);
        for(i = 0; i < num_groups; ++i){
            colours[i] = 0;
        }
        colours[colour] = 1;
    }
    std::cout << "Rank = " << world_rank << " : colour = " << colour << std::endl;
    MPI_Comm_split(MPI_COMM_WORLD, colour, world_rank, &comm);

    /* reader code */
    if(colour != num_groups){
	    comm_writer.push_back(MPI_COMM_NULL);
        for(i = 0; i < colour; ++i){
            MPI_Comm_split(MPI_COMM_WORLD, colours[i], world_rank, &comm_writer[0]);
        }
        MPI_Comm_split(MPI_COMM_WORLD, colours[colour], world_rank, &comm_reader);
        enginePair = "globalArray_";
	    enginePair.append(std::to_string(colour));
	    std::replace(enginePair.begin(), enginePair.end(), '/', '_');
        std::cout << "Local rank: " << world_rank << ": Reader "<< enginePair <<std::endl;
        adios_catalyst_init(comm, comm_reader, enginePair, 0);
        for(i = colour+1; i < num_groups; ++i){
            MPI_Comm_split(MPI_COMM_WORLD, colours[i], world_rank, &comm_writer[0]);
        }
        adios_catalyst();
        /*After first step of the first writer and reader pair, rank 0 in reader would send back to writer 0 to inform it this reader can take new job again.*/
    }else{ 
    /*Here starts the writer code*/
        comm_f = MPI_Comm_c2f(comm);
	    nek_init_malleable_insitu_(&comm_f);
        for(i = 0; i < num_groups; ++i){
            nek_solve_malleable_insitu_first_(&i);
	        comm_writer.push_back(MPI_COMM_NULL);
            MPI_Comm_split(MPI_COMM_WORLD, colours[i], world_rank, &comm_writer[i]);
            enginePair = "globalArray_";
            enginePair.append(std::to_string(i));
            std::replace(enginePair.begin(), enginePair.end(), '/', '_');
            std::cout << "Local rank: " << world_rank << ": Writer "<< enginePair << std::endl;
	    int worldComm_f = MPI_Comm_c2f(comm_writer[i]);
		   
            adios_writer_init(comm, comm_writer[i], enginePair);
            in_situ_init_(&worldComm_f);
	}
        nek_solve_malleable_insitu_(&num_groups);
		nek_end_();
    }

    return 0;
    
}

