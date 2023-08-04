#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include "3rd_party/adios_catalyst_malleable_insitu.h"
#include "3rd_party/nek_adios2.h"
extern "C"{
	void nek_init_malleable_insitu_(int * comm);
	void nek_solve_malleable_insitu_first_();
	void in_situ_init_(int * worldComm);
	void in_situ_check_();
	void nek_solve_malleable_insitu_();
	void nek_end_();
}


void free_string_array(char **array, int size){
    for(int i = 0; i < size; i++){
        free(array[i]);
    }
    free(array);
}


int main(int argc, char* argv[]){
	MPI_Group group = MPI_GROUP_NULL;
	MPI_Session session = MPI_SESSION_NULL;
	MPI_Comm comm = MPI_COMM_NULL;
	MPI_Comm comm_ = MPI_COMM_NULL;
	MPI_Comm worldComm = MPI_COMM_NULL;
	MPI_Info info = MPI_INFO_NULL;
	char main_pset[MPI_MAX_PSET_NAME_LEN];
	char boolean_string[16], nprocs[] = "2", **input_psets, **output_psets, host[64];  
	int original_rank, new_rank, flag = 0, dynamic_process = 0, noutput, op;
	int comm_f, worldComm_f;
	int world_rank, world_size;

	gethostname(host, 64);
	char *dict_key = strdup("main_pset"); // The key used to store the name of the new main PSet in the PSet Dictionary

	/* We start with the mpi://WORLD PSet as main PSet */
	strcpy(main_pset, "mpi://WORLD");

	/* Initialize the MPI Session */
	MPI_Session_init(MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &session);

	/* Get the info from our mpi://WORLD pset */
	MPI_Session_get_pset_info (session, main_pset, &info);

	/* get value for the 'mpi_dyn' key -> if true, this process was added dynamically */
	MPI_Info_get(info, "mpi_dyn", 6, boolean_string, &flag);
	MPI_Info_free(&info);

	MPI_Group_from_session_pset (session, main_pset, &group);
        MPI_Comm_create_from_group(group, "mpi.forum.example", MPI_INFO_NULL, MPI_ERRORS_RETURN, &comm);
        printf("Wrong rank %d!\n", original_rank);
	MPI_Comm_rank(comm, &original_rank);
	printf("Rank %d!\n", original_rank);
        MPI_Group_free(&group);

	/* if mpi://WORLD is a dynamic PSet retrieve the name of the main PSet stored on mpi://WORLD */
	if(dynamic_process = (flag && 0 == strcmp(boolean_string, "True"))){
		/* Lookup the value for the "main_pset" key in the PSet Dictionary and use it as our main PSet */
		MPI_Session_get_pset_data (session, main_pset, main_pset, (char **) &dict_key, 1, true, &info);
		MPI_Info_get(info, "main_pset", MPI_MAX_PSET_NAME_LEN, main_pset, &flag);
		MPI_Info_free(&info);
	}


	/* create a communcator from our main PSet */
	MPI_Group_from_session_pset (session, main_pset, &group);
	MPI_Comm_create_from_group(group, "mpi.forum.example", MPI_INFO_NULL, MPI_ERRORS_RETURN, &worldComm);
	MPI_Comm_rank(worldComm, &original_rank);
	MPI_Group_free(&group);
	printf("Rank %d: Host = '%s', Main PSet = '%s'. I am '%s'!\n", original_rank, host, main_pset, dynamic_process ? "dynamic" : "original");

	/* Original processes will switch to a grown communicator */
	if(!dynamic_process){
		/* One process needs to request the set operation and publish the kickof information */
		//MPI_Comm_dup(worldComm, &comm);
		comm_f = MPI_Comm_c2f(comm);
		nek_init_malleable_insitu_(&comm_f);
		nek_solve_malleable_insitu_first_();
		if(original_rank == 0){
			
			/* Request the GROW operation */
			op = MPI_PSETOP_GROW;
			
			/* We add nprocs = 2 processes*/
			MPI_Info_create(&info);
			MPI_Info_set(info, "mpi_num_procs_add", nprocs);
			
			/* Thein PSet is the input PSet of the operation */
			input_psets = (char **) malloc(1 * sizeof(char*));
			input_psets[0] = strdup(main_pset);
			noutput = 0;
			
			/* Send the Set Operation request */
			MPI_Session_dyn_v2a_psetop(session, &op, input_psets, 1, &output_psets, &noutput, info);
			MPI_Info_free(&info);
			printf("OP=%d\n", op);
			printf("4%s\n", output_psets[0]);
			printf("4%s\n", output_psets[1]);
			/* Publish the name of the new main PSet on the delta Pset */
			MPI_Info_create(&info);
			printf("40\n");
			MPI_Info_set(info, "main_pset", output_psets[1]);
			printf("41\n");
			MPI_Session_set_pset_data(session, output_psets[0], info);
			printf("42\n");
			MPI_Info_free(&info);
			printf("43\n");
			free_string_array(input_psets, 1);
			printf("44\n");
			free_string_array(output_psets, noutput);
			printf("5\n");
		}

		/* All processes can query the information about the pending Set operation */      
		MPI_Session_dyn_v2a_query_psetop(session, main_pset, main_pset, &op, &output_psets, &noutput);

		/* Lookup the name of the new main PSet stored on the delta PSet */
		MPI_Session_get_pset_data (session, main_pset, output_psets[0], (char **) &dict_key, 1, true, &info);
		MPI_Info_get(info, "main_pset", MPI_MAX_PSET_NAME_LEN, main_pset, &flag); 
		free_string_array(output_psets, noutput);
		MPI_Info_free(&info);

		/* Disconnect from the old communicator */
    		MPI_Comm_disconnect(&worldComm);
		
		/* create a new ommunicator from the new main PSet*/
		MPI_Group_from_session_pset (session, main_pset, &group);
		MPI_Comm_create_from_group(group, "mpi.forum.example", MPI_INFO_NULL, MPI_ERRORS_RETURN, &worldComm);
		MPI_Comm_rank(worldComm, &new_rank);
		MPI_Group_free(&group);
		
		/* Indicate completion of the Pset operation*/
		if(original_rank == 0){
			MPI_Session_dyn_finalize_psetop(session, main_pset);
		}

		printf("Rank %d: Host = '%s', Main PSet = '%s'. I am 'original'!\n", new_rank, host, main_pset);
		worldComm_f = MPI_Comm_c2f(worldComm);
    		MPI_Comm_rank(worldComm, &world_rank);
    		MPI_Comm_size(worldComm, &world_size);
    		MPI_Comm_split(worldComm, 1024, world_rank, &comm_);
		
		adios_MPI_set(comm, worldComm);
		in_situ_init_(&worldComm_f);
		in_situ_check_();
		nek_solve_malleable_insitu_();
		nek_end_();
	}else{
		MPI_Comm_rank(worldComm, &world_rank);
    		MPI_Comm_size(worldComm, &world_size);
    		MPI_Comm_split(worldComm, 2048, world_rank, &comm_);
		adios_catalyst(comm, worldComm);
	}


    MPI_Barrier(comm);
    /* Disconnect from the old communicator */
    MPI_Comm_disconnect(&comm);

    MPI_Barrier(worldComm);
    MPI_Comm_disconnect(&worldComm);

    /* Finalize the MPI Session */
    MPI_Session_finalize(&session);

    return 0;

}



