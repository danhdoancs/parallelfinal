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
void backSubstitutionPhase();
void debug_array(double* arr, int size, string msg);
void debug(string msg, double varr);

//vars
long str_len;
int noOfProcesses;
int myId;
double s_time, e_time;

const int N = 8;
const int M = N - 1;
const int n = log2(N);

double* A = new double[M] {
    0, 5, 1, 3, 6, 7, 2
};
double* B = new double[M] {
    1, 2, 9, 3, 5, 8, 4
};
double* C = new double[M] {
    3, 5, 1, 3, 4, 2, 0
};
double* D = new double[M] {
    6, 2, 4, 1, 5, 7, 8
};
double* x = new double[M];
double** P = new double*[4];
double* Pnext;
double* Pprevious;

int h;

/*
 * 
 */
int main(int argc, char** argv) {

    //communication
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &noOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

    //check if noOfProcess is not equal n
    if (noOfProcesses != N) {
        cout << "Please use " << N << "processors\n";
        return -1;
    }

    // wait for all to come together
    MPI_Barrier(MPI_COMM_WORLD);
    s_time = MPI_Wtime();

    //init P0
    P[0] = new double[4] {
        A[myId], B[myId], C[myId], D[myId]
    };

    //run
    reductionPhase();

    backSubstitutionPhase();

    MPI_Barrier(MPI_COMM_WORLD);
    e_time = MPI_Wtime();

    if (myId != 0)
        printf("Test: process %d received message %ld bytes in %f seconds\n", myId, str_len, e_time - s_time);

    MPI_Finalize();

    //debug_array(x, M, "array X");
    return 0;
}

void backSubstitutionPhase() {
    //    int k = pow(2, n - 1);
    //    x[k] = D[n - 1][k] / B[n - 1][k];
    //
    //    for (int k = n - 1; k > 0; k--) {
    //        int h = pow(2, k - 1);
    //        for (int i = h; i < pow(2, n); i += 2 * h) {
    //            if (myId == i)
    //                x[i] = (D[k - 1][i] - A[k - 1][i] * x[i - h] - C[k - 1][i] * x[i + h]) / (B[k - 1][i]);
    //        }
    //    }
}

void reductionPhase() {
    for (int j = 1; j < n; j++) {

        h = pow(2, j - 1);
        Pnext = new double[4];
        Pprevious = new double[4];

        if (myId - h >= 0) {

            //send values to other nodes
        }

        //get values from other nodes
        if (myId + h <= M) {

        }
        for (int i = pow(2, j); i < pow(2, n); i += pow(2, j)) {

            //debug("Reduction: i", (double)i)            
            if (myId == i - h) {
                printf("Processor %i send to %i\n", myId, i);
                MPI_Send(P[j - 1], 4, MPI_DOUBLE, myId, 0, MPI_COMM_WORLD);
            } else if (myId == i + h) {
                printf("Processor %i send to %i\n", myId, i);
                MPI_Send(P[j - 1], 4, MPI_DOUBLE, myId + h, 0, MPI_COMM_WORLD);
            } else if (myId == i) {

                printf("Processor %i receive from %i and %i\n", myId, myId - h, myId + h);
                MPI_Recv(Pprevious, 4, MPI_DOUBLE, myId - h, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(Pnext, 4, MPI_DOUBLE, myId + h, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //                debug("Reduction: myId", myId);
                computePi(j, i);

                //send values to other nodes
                //MPI_Send(greeting, strlen(greeting) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

                break;
            }
        }
    }
}

void computePi(int j, int i) {
    //init
    double e, f;

    //debug("(B[j - 1][i - h])", (B[j - 1][i - h]));
    //debug_array(B[j - 1], M, "B");

    printf("Processor %i computePi(j:%i, i:%i, h:%i, i-h:%i, i+h:%i)\n", i, j, i, h, i - h, i + h);
    debug_array(P[j - 1], 4, "P " + std::to_string(j - 1));

    if (i <= 0 or i >= N) {
        P[j] = new double[4] {
            0, 1, 0, 0
        };
        x[i] = 0;
    } else {



        debug_array(Pprevious, 4, "Pprevious ");
        debug_array(Pnext, 4, "Pnext ");

        e = -(P[j - 1][0]) / (Pprevious[1]);
        f = -(P[j - 1][2]) / (Pnext[1]);
        P[j][0] = e * Pprevious[0];
        P[j][2] = f * Pnext[2];
        P[j][1] = P[j - 1][1] + e * Pprevious[2] + f * Pnext[0];
        P[j][3] = P[j - 1][3] + e * Pprevious[3] + f * Pnext[3];
    }


    debug_array(P[j], M, "P " + std::to_string(j));

    printf("Processor %i computePi(j:%i, i:%i, h:%i, e:%f)\n", i, j, i, h, e);
}

//debug functions

void debug(const char *msg, double varr) {

    cout << "\nDebugging....." << msg << ": <" << varr << ">\n";
    fflush(stdout);
}

//debug a array

void debug_array(double* arr, int size, string msg = "") {
    int i = 0;
    cout << "\nDebugging..." << msg << ": ";

    for (; i < size; i++)
        cout << arr[i] << "  ";
    cout << "\n";

    fflush(stdout);
}

void array_show(int* arr, int size) {
    for (int i = 0; i < size; i++)
        cout << arr[i] << "  ";
}