/* 
 * File:   main.cpp
 * Author: david
 *
 * Created on April 15, 2015, 10:05 PM
 */

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>
#include <string.h> 
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;
void reductionPhase();
void computePi(int j, int i);

//vars
long str_len;
int noOfProcesses;
int myId;
double s_time, e_time;

const int N = 8;
const int M = N - 1;
const int n = log(N);

double** A = new double*[n];
double** B = new double*[n];
double** C = new double*[n];
double** D = new double*[n];
double* x = new double[M];

/*
 * 
 */
int main(int argc, char** argv) {

    //init first vars
    A[0] = new double[M] {
        0, 5, 1, 3, 6, 7, 2
    };

    B[0] = new double[M] {
        1, 2, 9, 3, 5, 8, 4
    };

    C[0] = new double[M] {
        3, 5, 1, 3, 4, 2, 0
    };
    D[0] = new double[M] {
        6, 2, 4, 1, 5, 7, 8
    };

    //communication
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &noOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

    //check if noOfProcess is not equal n
    if(noOfProcesses != M)
    {
        cout << "Please use " << M << "processors\n";
        return;
    }

    // wait for all to come together
    MPI_Barrier(MPI_COMM_WORLD);
    s_time = MPI_Wtime();

    //run
    reductionPhase();

    MPI_Barrier(MPI_COMM_WORLD);
    e_time = MPI_Wtime();

    if (myId != 0)
        printf("Test: process %d received message %ld bytes in %f seconds\n", myId, str_len, e_time - s_time);

    MPI_Finalize();
    return 0;
}

void backSubstitutionPhase() {
    int k = pow(2,n-1);
    x[k] = D[n-1][k] / B[n-1][k];
    
    for(int k = n-1; k > 0; k--)
    {
        int h = pow(2,k-1);
        for(int i = h; i < pow(2,n); i+= 2*h)
        {
            x[i] = 0;
        }
    }
}

void reductionPhase() {
    for (int j = 1; j < n; j++) {
        //init
        A[j] = new double[M];
        B[j] = new double[M];
        C[j] = new double[M];
        D[j] = new double[M];

        for (int i = pow(2, j); i < pow(2, n); i += pow(2, j)) {
            //processor 0 inits
            if (myId == 0) {              

            } else {

                //MPI_Recv(msg, str_len, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
            if(myId == i)
                computePi(j, i);
        }
    }
}

void computePi(int j, int i) {
    //init
    int h = pow(2, j - 1);

    if (i <= 0 or i >= N) {
        A[j][i] = 0;
        B[j][i] = 1;
        C[j][i] = 0;
        D[j][i] = 0;
        x[i] = 0;
    } else {
        double e = -(A[j - 1][i]) / (B[j - 1][i - h]);
        double f = -(C[j - 1][i]) / (B[j - 1][i - h]);
        A[j][i] = e * A[j - 1][i - h];
        C[j][i] = f * C[j - 1][i + h];
        B[j][i] = B[j - 1][i] + e * C[j - 1][i - h] + f * A[j - 1][i + h];
        D[j][i] = D[j - 1][i] + e * D[j - 1][i - h] + f * D[j - 1][i + h];
    }
}

