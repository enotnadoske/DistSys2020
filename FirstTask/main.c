#include <mpi.h>
#include <stdio.h>
#include "mpi-ext.h"

bool shift1 = false;

void modify_DSM(int * a,int rank, int max_rank, int coord){
    //rank == coord rank
    if(rank == coord){
        MPI_Status st;
        int value[3];
        value[0] = 1;
        value[1] = rank;
        value[2] = 0;//modify operation
        for(int i = 0; i < max_rank; i++){
            if (i == rank) continue;
            MPI_Send(&value, 3, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        a[rank] = 1;
        //printf("shiiir");
    }else{
        printf("shiiit");
        int value[3];
        value[0] = 1;
        value[1] = rank;
        value[2] = 1;//coordinate operation
        MPI_Send(&value, 3, MPI_INT, coord, 0, MPI_COMM_WORLD);
        //a[rank] = 1;
    }    
}


int main(int argc, char** argv) {
    //array models DSM
    int a[10] = {0,0,0,0,0,0,0,0,0,0}; 
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    if (argc == 2){
        //BAD SITUATION
        shift1 = true;// lets assume that coordinators are shifted to the right by 1 
                      // so to the first cell corresponds second process as coordinator 
    }

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Status st;
    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    
    if (!shift1){
        modify_DSM(a,world_rank, world_size,world_rank);
        for (int i = 0; i < world_size - 1; i++ ){
            int value[3];
            MPI_Recv(&value, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&st);
            a[value[1]] = value[0];
        }
    }else{
        if (world_rank == 9){
            modify_DSM(a,world_rank, world_size,0);
        }
        else{
            modify_DSM(a,world_rank, world_size,world_rank+1);
        }
        for (int i = 0; i < world_size; i++ ){
            int value[3];
            MPI_Recv(&value, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&st);
            if(value[2] == 0){
                a[value[1]] = value[0];
            }
            else{
                int tmp[3];
                tmp[0] = value[0];
                tmp[1] = value[1];
                tmp[2] = 0;
                //printf("%d\n", a[2]);
                //fflush(stdout);
                for(int i = 0; i < world_size; i++){
                    if (i == world_rank) continue;
                    MPI_Send(&tmp, 3, MPI_INT, i, 0, MPI_COMM_WORLD);
                }
                a[value[1]] = value[0];
            }
        }
    }        

    // Print off a hello world message
   // printf("Hello world from processor %s, rank %d out of %d processors\n",
     //      processor_name, world_rank, world_size);
    for(int i = 0; i < 10 ; i++){
        printf("%d",a[i]);
    }
    printf("\n");

    // Finalize the MPI environment.
    MPI_Finalize();
}
