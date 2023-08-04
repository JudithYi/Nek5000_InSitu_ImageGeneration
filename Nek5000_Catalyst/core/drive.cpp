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
    MPI_Comm comm_reader = MPI_COMM_NULL;
    std::vector<MPI_Comm> comm_writer;
    MPI_Info info = MPI_INFO_NULL;
    char main_pset[MPI_MAX_PSET_NAME_LEN];
    char delta_pset[MPI_MAX_PSET_NAME_LEN];
    char old_main_pset[MPI_MAX_PSET_NAME_LEN];
    char _keys[][MPI_MAX_INFO_KEY] = {"next_main_pset"}; // Name of the main PSet to use
    char *keys[1] = {_keys[0]};
    char boolean_string[16], **input_psets, **output_psets, host[64];
    char nprocs[] = "36";  
    int original_rank, original_size, new_size, flag = 0, dynamic_process = 0, noutput, op;
    int flag_loop = 0;
    int i = 0;
    int rc;
    int firstPair = 0;
    MPI_Request request;
    MPI_Status MPI_status;
    std::string enginePair;
    int comm_f;

    int readerDone = 0;

    int n = 0;

    /* Parse the command line arguments*/
    gethostname(host, 64);

    /* We start with the mpi://WORLD PSet as main PSet */
    strcpy(main_pset, "mpi://WORLD");
    strcpy(delta_pset, "");
    strcpy(old_main_pset, "mpi://WORLD");

    /* Initialize the MPI Session */
    MPI_Session_init(MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &session);

    /* Get the info from our mpi://WORLD pset */
    MPI_Session_get_pset_info (session, main_pset, &info);

    /* get value for the 'mpi_dyn' key -> if true, this process was added dynamically */
    MPI_Info_get(info, "mpi_dyn", 6, boolean_string, &flag);
    MPI_Info_free(&info);

    MPI_Group_from_session_pset (session, main_pset, &group);
    MPI_Comm_create_from_group(group, "mpi.forum.example", MPI_INFO_NULL, MPI_ERRORS_RETURN, &comm);
    MPI_Comm_rank(comm, &original_rank);
    MPI_Comm_size(comm, &original_size);
    MPI_Group_free(&group);

    /* if mpi://WORLD is a dynamic PSet retrieve the name of the main PSet stored on mpi://WORLD */
    if(dynamic_process = (flag && 0 == strcmp(boolean_string, "True"))){
        /* Lookup the value for the "main_pset" key in the PSet Dictionary and use it as our main PSet */
        MPI_Session_get_pset_data (session, main_pset, main_pset, (char **) &keys[0], 1, true, &info);
        MPI_Info_get(info, "next_main_pset", MPI_MAX_PSET_NAME_LEN, main_pset, &flag);
        if(!flag){
            printf("No 'next_main_pset' was provided for dynamic process. Terminate.\n");
            MPI_Session_finalize(&session);
            return -1;
        }
        MPI_Info_free(&info);

        /* create communication communcator (comm_reader) from our main PSet */
        MPI_Group_from_session_pset (session, main_pset, &group);
	    rc=MPI_Comm_create_from_group(group, "mpi.forum.examples", MPI_INFO_NULL, MPI_ERRORS_RETURN, &comm_reader);
	    if(rc!=MPI_SUCCESS){
		    printf("MPI_comm_create_from_group failed with rc = %d\n", rc);
        	return rc;
	    }
	    MPI_Comm_size(comm_reader, &new_size);
        MPI_Group_free(&group);
        /* start_ins-situ(in-situ_comm, communication_comm); */
        if(!original_rank) printf("Host = '%s', Main PSet = '%s'. I am '%s'! Local communicator size is: %d, communication communicator size is %d. \n", host, main_pset, dynamic_process ? "dynamic" : "original", original_size, new_size);
        /*Check if this is the first pair of writer and reader*/
        if(!original_rank) MPI_Recv(&firstPair, 1, MPI_INT, 0, 0, comm_reader, &MPI_status);

        enginePair = "globalArray_";
	    enginePair.append(main_pset);
	    std::replace(enginePair.begin(), enginePair.end(), '/', '_');
        if(!original_rank) std::cout << "Local rank: " << original_rank << ": Reader "<< enginePair <<std::endl;

        /*After first step of the first writer and reader pair, rank 0 in reader would send back to writer 0 to inform it this reader can take new job again.*/
        adios_catalyst(comm, comm_reader, enginePair, firstPair); 
	/*
        readerDone = 1;
        if(firstPair){
            if(!original_rank) MPI_Send(&readerDone, 1, MPI_INT, 0, 0, comm_reader);
        }
        adios_catalyst_run();
	*/
    }else{ 
    /*Here starts the writer code*/
        comm_f = MPI_Comm_c2f(comm);
	nek_init_malleable_insitu_(&comm_f);
        nek_solve_malleable_insitu_first_(&i);
	    /*start_sim(sim_comm)*/
        i=0;
        /*The first pair of wrtier and reader*/
        if(original_rank == 0){
            /* Request the GROW operation */
            op = MPI_PSETOP_GROW;
            
            /* We add nprocs = 2 processes*/
            MPI_Info_create(&info);
            MPI_Info_set(info, "mpi_num_procs_add", nprocs);
            
            /* The main PSet is the input PSet of the operation */
            input_psets = (char **) malloc(1 * sizeof(char*));
            input_psets[0] = strdup(old_main_pset);
            noutput = 0;
            
            /* Send the Set Operation request */
            MPI_Session_dyn_v2a_psetop(session, &op, input_psets, 1, &output_psets, &noutput, info);
            MPI_Info_free(&info);
            
            /* Publish the name of the new main PSet on the delta Pset */
            if(MPI_PSETOP_NULL != op){

                MPI_Info_create(&info);
                MPI_Info_set(info, "next_main_pset", output_psets[1]);
                MPI_Session_set_pset_data(session, output_psets[0], info);
                MPI_Info_free(&info);
            }else{
                printf("op=%d.\n", op);
            }

            free_string_array(input_psets, 1);
            free_string_array(output_psets, noutput);

        }
        /* All processes can query the information about the pending Set operation */        
        MPI_Session_dyn_v2a_query_psetop(session, main_pset, old_main_pset, &op, &output_psets, &noutput);

        /* Lookup the name of the new main PSet stored on the delta PSet */
        MPI_Session_get_pset_data (session, main_pset, output_psets[0], (char **) &keys[0], 1, true, &info);
        MPI_Info_get(info, "next_main_pset", MPI_MAX_PSET_NAME_LEN, main_pset, &flag); 
        free_string_array(output_psets, noutput);
        if(!flag){
            printf("could not find next_main_pset on PSet %s. This should never happen! Terminate.\n", main_pset);
            return -1;
        }
        /* create a new ommunicator from the new main PSet*/
	comm_writer.push_back(MPI_COMM_NULL);
        MPI_Info_free(&info);
        MPI_Group_from_session_pset (session, main_pset, &group);
        rc=MPI_Comm_create_from_group(group, "mpi.forum.examples", MPI_INFO_NULL, MPI_ERRORS_RETURN, &comm_writer[i]);
        if(rc!=MPI_SUCCESS){
            printf("MPI_comm_create_from_group failed with rc = %d\n", rc);
            return rc;
        }
        MPI_Comm_size(comm_writer[i], &new_size);
        MPI_Group_free(&group);
        /* Indicate completion of the Pset operation*/
        if(original_rank == 0){
            MPI_Session_dyn_finalize_psetop(session, main_pset);
        }
        if(!original_rank) printf("%d: Host = '%s', Main PSet = '%s'. I am '%s'! Local communicator size is: %d, communication communicator size is %d. \n", i, host, main_pset, dynamic_process ? "dynamic" : "original", original_size, new_size);
        firstPair = 1;
        if(!original_rank) MPI_Send(&firstPair, 1, MPI_INT, original_size, 0, comm_writer[i]);
        enginePair = "globalArray_";
        enginePair.append(main_pset);
        std::replace(enginePair.begin(), enginePair.end(), '/', '_');
        if(!original_rank) std::cout << "Local rank: " << original_rank << ": Writer "<< enginePair << std::endl;
        
        adios_writer_init(comm, comm_writer[i], enginePair);
		int worldComm_f = MPI_Comm_c2f(comm_writer[i]);
		in_situ_init_(&worldComm_f);
        
        if(!original_rank) MPI_Irecv(&readerDone, 1, MPI_INT, original_size, 0, comm_writer[i], &request);
        if(!original_rank) MPI_Test(&request, &flag_loop, &MPI_status);
        ++i;
        ++n;
        while(!flag_loop){
            /*PSetOp(output_names);*/
            nek_solve_malleable_insitu_first_(&i);
            /* One process needs to request the set operation and publish the kickof information */
            if(original_rank == 0){

                /* Request the GROW operation */
                op = MPI_PSETOP_GROW;
                
                /* We add nprocs = 2 processes*/
                MPI_Info_create(&info);
                MPI_Info_set(info, "mpi_num_procs_add", nprocs);
                
                /* The main PSet is the input PSet of the operation */
                input_psets = (char **) malloc(1 * sizeof(char*));
                input_psets[0] = strdup(old_main_pset);
                noutput = 0;
                
                /* Send the Set Operation request */
                MPI_Session_dyn_v2a_psetop(session, &op, input_psets, 1, &output_psets, &noutput, info);
                MPI_Info_free(&info);
                
                /* Publish the name of the new main PSet on the delta Pset */
                if(MPI_PSETOP_NULL != op){

                    MPI_Info_create(&info);
                    MPI_Info_set(info, "next_main_pset", output_psets[1]);
                    MPI_Session_set_pset_data(session, output_psets[0], info);
                    MPI_Info_free(&info);
                }else{
                    printf("op=%d.\n", op);
                }

                free_string_array(input_psets, 1);
                free_string_array(output_psets, noutput);

            }
            /* All processes can query the information about the pending Set operation */        
            MPI_Session_dyn_v2a_query_psetop(session, main_pset, old_main_pset, &op, &output_psets, &noutput);

            /* Lookup the name of the new main PSet stored on the delta PSet */
            MPI_Session_get_pset_data (session, main_pset, output_psets[0], (char **) &keys[0], 1, true, &info);
            MPI_Info_get(info, "next_main_pset", MPI_MAX_PSET_NAME_LEN, main_pset, &flag); 
            free_string_array(output_psets, noutput);
            if(!flag){
                printf("could not find next_main_pset on PSet %s. This should never happen! Terminate.\n", main_pset);
                return -1;
            }
            /* create a new ommunicator from the new main PSet*/
            MPI_Info_free(&info);
            MPI_Group_from_session_pset (session, main_pset, &group);
            rc=MPI_Comm_create_from_group(group, "mpi.forum.examples", MPI_INFO_NULL, MPI_ERRORS_RETURN, &comm_writer[i]);
            if(rc!=MPI_SUCCESS){
                printf("MPI_comm_create_from_group failed with rc = %d\n", rc);
                return rc;
            }
	        MPI_Comm_size(comm_writer[i], &new_size);
	        MPI_Group_free(&group);
            /* Indicate completion of the Pset operation*/
            if(original_rank == 0){
                MPI_Session_dyn_finalize_psetop(session, main_pset);
            }
            printf("%d: Host = '%s', Main PSet = '%s'. I am '%s'! Local communicator size is: %d, communication communicator size is %d. \n", i, host, main_pset, dynamic_process ? "dynamic" : "original", original_size, new_size);

            firstPair = 0;
            if(!original_rank) MPI_Send(&firstPair, 1, MPI_INT, original_size, 0, comm_writer[i]);
            enginePair = "globalArray_";
	        enginePair.append(main_pset);
	        std::replace(enginePair.begin(), enginePair.end(), '/', '_');
            std::cout << "Local rank: " << original_rank << ": Writer "<< enginePair << std::endl;
            adios_writer_init(comm, comm_writer[i], enginePair);
		    int worldComm_f = MPI_Comm_c2f(comm_writer[i]);
		    in_situ_init_(&worldComm_f);

            ++i;
            ++n;
            if(!original_rank) MPI_Test(&request, &flag_loop, &MPI_status);
            MPI_Bcast(&flag_loop, 1, MPI_INT, 0, comm);
        }
	std::cout << original_rank << ": No more processors added n = " << n << std::endl;
        nek_solve_malleable_insitu_(&n);
		nek_end_();
    }
    MPI_Barrier(comm);

    /* Disconnect from the old communicator */
    MPI_Comm_disconnect(&comm);
    if(dynamic_process){
        MPI_Barrier(comm_reader);

        /* Disconnect from the old communicator */
        MPI_Comm_disconnect(&comm_reader);
    }else{
        for(i=0;i<n;++i){
            MPI_Barrier(comm_writer[i]);

            /* Disconnect from the old communicator */
            MPI_Comm_disconnect(&comm_writer[i]);
        }
    }

    /* Finalize the MPI Session */
    MPI_Session_finalize(&session);

    return 0;
    
}

