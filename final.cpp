//Parallel Programming
//Homework 2
//Author: David Doan

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>
#include <string.h>

using namespace std;

int main(void) {
    
    //vars
    long str_len;
    int noOfProcesses;
    int myId;
    double s_time, e_time;

    //communication
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &noOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

    //test different data size
    int no_k = 0;
    for (int k = 0; k <= 6; k = k + 2) {

        str_len = pow(10, k); // 1 char = 1 byte
        char msg[str_len];

        //test 20 times for each size
        for (int i = 1; i <= 20; i++) {

            // wait for all to come together
            MPI_Barrier(MPI_COMM_WORLD);
            s_time = MPI_Wtime();

            if (myId == 0) {

                sprintf(msg, "");

                for (int q = 1; q < noOfProcesses; q++) {
                    MPI_Send(msg, str_len, MPI_CHAR, q, 0, MPI_COMM_WORLD);
                }

            } else {

                MPI_Recv(msg, str_len, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            MPI_Barrier(MPI_COMM_WORLD);
            e_time = MPI_Wtime();

            if (myId != 0)
                printf("Test %d: process %d received message %ld bytes in %f seconds\n", no_k * 20 + i, myId, str_len, e_time - s_time);
        }

        no_k++;
    }

    MPI_Finalize();
    return 0;
}