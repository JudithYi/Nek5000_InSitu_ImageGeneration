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
	void nek_init_insitu_(int * comm);
	void in_situ_init_(int * type_index);
	void nek_solve_();
	void nek_end_();
}

int main(int argc, char* argv[]){

    MPI_Comm comm = MPI_COMM_NULL;
    MPI_Comm comm_reader = MPI_COMM_NULL;
    std::vector<MPI_Comm> comm_writer;
    int i, j;
    
    std::string enginePair;
    int colour;
    std::vector<int> group_start;
    std::vector<int> colours;
    int comm_f;
    int world_rank, world_size;
    int cores_per_socket = std::stoi(argv[1]);
    int type_in_situ = std::stoi(argv[2]);
    std::vector<int> num_groups(type_in_situ);
    std::vector<int> cores_in_situ(type_in_situ);
    std::vector<int> fre_in_situ(type_in_situ);
    int total_group = 0;
    int total_insitu_cores = 0;
    for(i = 0; i < type_in_situ; ++i){
        num_groups[i] = std::stoi(argv[3 + i * 3]);
        cores_in_situ[i] = std::stoi(argv[4 + i * 3]);
        fre_in_situ[i] = std::stoi(argv[5 + i * 3]);
        total_group += num_groups[i];
        total_insitu_cores += num_groups[i] * cores_in_situ[i];
    }
    colours.resize(total_group);
    group_start.resize(total_group+1);
    group_start[0] = 0;
    int tmp_group = 0;
    for(i = 0; i < type_in_situ; ++i){
        for(j = 0; j < num_groups[i]; ++j){
            group_start[tmp_group + j + 1] = group_start[tmp_group + j] + cores_in_situ[i];
        }
        tmp_group += num_groups[i];
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if(world_rank%cores_per_socket>=total_insitu_cores){
        colour = total_group;
        for(i = 0; i < total_group; ++i){
            colours[i] = 1;
        }
    }else{
        //colour = static_cast<int>(world_rank%cores_per_socket/cores_in_situ);
        for(i = 0; i < total_group; ++i){
            colours[i] = 0;
            if(world_rank%cores_per_socket>=group_start[i] && world_rank%cores_per_socket<group_start[i+1]) colour = i;
        }
        colours[colour] = 1;
    }
    std::cout << "Rank = " << world_rank << " : colour = " << colour << std::endl;
    MPI_Comm_split(MPI_COMM_WORLD, colour, world_rank, &comm);

    /* reader code */
    if(colour != total_group){
	    comm_writer.push_back(MPI_COMM_NULL);
        for(i = 0; i < colour; ++i){
            MPI_Comm_split(MPI_COMM_WORLD, colours[i], world_rank, &comm_writer[0]);
        }
        MPI_Comm_split(MPI_COMM_WORLD, colours[colour], world_rank, &comm_reader);
        enginePair = "globalArray_";
	    enginePair.append(std::to_string(colour));
	    std::replace(enginePair.begin(), enginePair.end(), '/', '_');
        std::cout << "Local rank: " << world_rank << ": Reader "<< enginePair <<std::endl;
        adios_catalyst_init(comm, comm_reader, enginePair, 0, colour);
        for(i = colour+1; i < total_group; ++i){
            MPI_Comm_split(MPI_COMM_WORLD, colours[i], world_rank, &comm_writer[0]);
        }
        adios_catalyst();
        /*After first step of the first writer and reader pair, rank 0 in reader would send back to writer 0 to inform it this reader can take new job again.*/
    }else{ 
    /*Here starts the writer code*/
        comm_f = MPI_Comm_c2f(comm);
	    nek_init_insitu_(&comm_f);
        init_multiple_type(type_in_situ);
        int tmp_group = 0;
        int tmp_idx;
        for(i = 0; i < type_in_situ; ++i){
            for(j = 0; j < num_groups[i]; ++j){
                tmp_idx = tmp_group + j;
                comm_writer.push_back(MPI_COMM_NULL);
                MPI_Comm_split(MPI_COMM_WORLD, colours[tmp_idx], world_rank, &comm_writer[tmp_idx]);
                enginePair = "globalArray_";
                enginePair.append(std::to_string(tmp_idx));
                std::replace(enginePair.begin(), enginePair.end(), '/', '_');
                std::cout << "Local rank: " << world_rank << ": Writer "<< enginePair << std::endl;
                adios_writer_init(comm, comm_writer[tmp_idx], enginePair, fre_in_situ[i], i);
                in_situ_init_(&i);
            }
            tmp_group += num_groups[i];
        }
        nek_solve_();
		nek_end_();
    }

    return 0;
    
}

