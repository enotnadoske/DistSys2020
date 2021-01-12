#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <stdio.h>
#include <sys/time.h>
#include <mpi.h>
#include "mpi-ext.h"
#include <signal.h>


void prt1a(char *t1, double *v, int n, char *t2);

void print_matrix(double *a);


MPI_Comm main_comm = MPI_COMM_WORLD;
int N ;
double *A;
#define A(i, j) A[(i) * (N + 1) + (j)]
double *X;
bool reverse_sub = false;
bool err_fl = false;
int proc_num, myrank;
//int errs = 0;
//bool was = false;

int additional_procs = 2;

//char** v;

//bool trouble_barrier = false;


static void saveData()
{
    if (myrank == 0) {
        FILE* f = fopen("elim_save.bin", "wb");
        fwrite(&A[0], sizeof(double),  (N-1)*N, f);
        fclose(f);

        if (reverse_sub) {
            FILE* f = fopen("rs_save.bin", "wb");
            fwrite(&X[0], sizeof(double),  N, f);
            fclose(f);
        }
    }
    //printf("here-\n");
    //fflush(stdout);
    //MPI_Barrier(main_comm);
    MPI_Barrier(main_comm);
}

static void loadData()
{

    
    //if(myrank < proc_num){
        FILE* f = fopen("elim_save.bin", "rb");
        fread(&A[0], sizeof(double), (N-1)*N, f);
        fclose(f);
        if (reverse_sub) {
            FILE* f = fopen("rs_save.bin", "wb");
            fwrite(&X[0], sizeof(double),  N, f);
            fclose(f);
        }
    //}    
    //printf("Proc %d\n", myrank);
    //fflush(stdout);
    MPI_Barrier(main_comm);
}


static void verbose_errhandler(MPI_Comm* pcomm, int* perr, ...)
{
    //errs++;
    MPI_Comm comm = *pcomm;
    int err = *perr;
    char errstr[MPI_MAX_ERROR_STRING];
    int i, rank, size, nf, len, eclass;
    MPI_Group group_c, group_f;
    int *ranks_gc, *ranks_gf;


    MPI_Error_class(err, &eclass);
    if( MPIX_ERR_PROC_FAILED != eclass ) {
        printf("Rank %d / %d: Notified of error %s. %d found dead: { ",
           rank, size, errstr, nf);
        MPI_Abort(comm, err);// stops all procesees if some other err occured
    }
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    //was = true;
    MPIX_Comm_failure_ack(comm);
    MPIX_Comm_failure_get_acked(comm, &group_f);
    MPI_Group_size(group_f, &nf);
    MPI_Error_string(err, errstr, &len);
    // MPIX_Comm_revoke(comm);//some error occures if revoke is used


    printf("Rank %d / %d: Notified of error %s. %d found dead: { ",
           rank, size, errstr, nf);

    ranks_gf = (int*)malloc(nf * sizeof(int));
    ranks_gc = (int*)malloc(nf * sizeof(int));
    MPI_Comm_group(comm, &group_c);
    for(i = 0; i < nf; i++)
        ranks_gf[i] = i;
    MPI_Group_translate_ranks(group_f, nf, ranks_gf,
                              group_c, ranks_gc);
    for(i = 0; i < nf; i++)
        printf("%d ", ranks_gc[i]);
    printf("}\n");
    additional_procs -= nf;
    //MPIX_Comm_revoke(main_comm);
    MPIX_Comm_shrink(comm, &main_comm);// leave out of communicator dead processes
    MPI_Comm_rank(main_comm, &myrank);
    MPI_Comm_size(main_comm, &proc_num);
    
    proc_num -= additional_procs;
    loadData();
    //fflush(stdout);
    err_fl = true;
    free(ranks_gc);
    free(ranks_gf);
    //printf("%d\n", myrank);
    //fflush(stdout);
    MPI_Barrier(main_comm);


}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    //* adding extra processes to calculate on the same number of processes if error occures*
    
    if(argc == 3){ //enable custum number of additional proc
        additional_procs = atoi(argv[2]);
    }
    MPI_Comm parent_comm, intercom_comm;
    MPI_Comm_get_parent(&parent_comm);
    if (parent_comm == MPI_COMM_NULL) {
        MPI_Comm_spawn("./a.out", 
                argv + 1, 
                additional_procs, 
                MPI_INFO_NULL, 0, 
                main_comm, 
                &intercom_comm, 
                MPI_ERRCODES_IGNORE);
    } else {
        intercom_comm = parent_comm;
        //printf("helll%d", atoi(argv[1]));
        fflush(stdout);

    }
    MPI_Intercomm_merge(intercom_comm, (parent_comm == MPI_COMM_NULL ? 0 : 1), &main_comm);

    //*intializtion of some variables*
    MPI_Comm_size(main_comm, &proc_num);
    MPI_Comm_rank(main_comm, &myrank);
    proc_num = proc_num - additional_procs;
    fflush(stdout);
    MPI_Errhandler errh;
    MPI_Comm_create_errhandler(verbose_errhandler, &errh);
    MPI_Comm_set_errhandler(main_comm, errh);//error handler setup completed
    int i, j, k;
    bool flag = false;
    N = atoi(argv[1]);//matrix size

    fflush(stdout);
    /* create arrays */
    A = (double *)malloc(N * (N + 1) * sizeof(double));
    X = (double *)malloc(N * sizeof(double));
    if (myrank == 0) printf("GAUSS %dx%d\n----------------------------------\n", N, N);
    srand(12345);
    /* initialize array A*/
    for (i = 0; i <= N - 1; i++)
        for (j = 0; j <= N; j++)
            if (i == j || j == N)
                A(i, j) = 1.f;
            else
                A(i, j) = 0.f;

    double time0 = MPI_Wtime();
    saveData();
    // if(!trouble_barrier){
    // }else{
    //     trouble_barrier = false;
    // }

    /* elimination */

    
    //err_fl = false;
    if (myrank == 0){
        MPI_Barrier(main_comm);
        raise(SIGKILL);
    }else{
        MPI_Barrier(main_comm);
    }
    //bool was = false; 
    flag = true;
    while (err_fl || flag) {
        err_fl = false;
        
            for (i = 0; i < N - 1; i++)
            {
                //printf("myrank = %d\n", myrank);
                MPI_Bcast(&A(i, i + 1), N - i, MPI_DOUBLE, i % proc_num, main_comm);
                if(myrank >= proc_num){ continue;}
                for (k = myrank; k <= N - 1; k += proc_num) {
                    if (k < i + 1) {
                        continue;
                    }
                    for (j = i + 1; j <= N; j++) {
                        A(k, j) = A(k, j) - A(k, i) * A(i, j) / A(i, i);
                    }
                }
            }
        flag = false;
        if (!err_fl)
                //printf(" %d HERE4!!!!1 , %d\n", myrank, err_fl);
            MPI_Barrier(main_comm);

    }
    //printf(" %d HERE4!!!!1 , %d\n", myrank, err_fl);

    saveData();
    //printf(" %d HERE4!!!!1 , %d\n", myrank, err_fl);

    
    /* reverse substitution */
    //printf("here2\n");
    //err_fl = false;
    flag = true;
    while (err_fl || flag) {
        err_fl = false;
        X[N - 1] = A(N - 1, N) / A(N - 1, N - 1);
        flag = false;
        if (!err_fl)
            MPI_Barrier(main_comm);
        //err_fl = false;
    }
    saveData();
    reverse_sub = true;

    if (myrank == 0){
        MPI_Barrier(main_comm);
        raise(SIGKILL);
    }else{
        MPI_Barrier(main_comm);
    }
    //err_fl = false;
    flag = true;
    while (err_fl || flag) {
        err_fl = false;  
            for (j = N - 2; j >= 0; j--) {
                if(myrank < proc_num){
                    for (k = myrank; k <= j; k += proc_num){
                        A(k, N) = A(k, N) - A(k, j + 1) * X[j + 1];
                    }
                }
                //printf("!!\n");    
                MPI_Bcast(&A(j, N), 1, MPI_DOUBLE, j % proc_num, main_comm);
                //printf("!!!!\n");
                X[j] = A(j, N) / A(j, j);
                if(err_fl){
                    break;
                }
            }
        //MPI_Barrier(main_comm);
        flag = false;
        if (!err_fl)
            MPI_Barrier(main_comm);
      
    }
    //printf("!!!!!");
    saveData();
    double time1 = MPI_Wtime();

    if (myrank == 0) {
        printf("Time in seconds=%gs\n", time1 - time0);
        prt1a("X=(", X, N > 9 ? 9 : N, "...)\n");
        printf("Final_proc_num - %d\n", proc_num );
    }

    free(A);
    free(X);
    
    MPI_Barrier(main_comm);
    MPI_Finalize();
    return 0;
}

void prt1a(char *t1, double *v, int n, char *t2)
{
    int j;
    printf("%s", t1);
    for (j = 0; j < n; j++)
        printf("%.4g%s", v[j], j % 10 == 9 ? "\n" : ", ");
    printf("%s", t2);
}

void print_matrix(double* a) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N + 1; j++)
            printf("%lf ", A(i, j));
        printf("\n");
    }
    printf("\n");
}


