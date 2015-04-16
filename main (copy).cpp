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

/*
 * 
 */
int main(int argc, char** argv) {

    //vars
    long str_len;
    int noOfProcesses;
    int myId;
    double s_time, e_time;
    int N;

    //communication
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &noOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);


    // wait for all to come together
    MPI_Barrier(MPI_COMM_WORLD);
    s_time = MPI_Wtime();

    str_len = 10;
    //char* msg;
    //processor 0 inits
    if (myId == 0) {

        //read file
        string line;
        ifstream myfile("input.txt");
        if (myfile.is_open()) {
            //get N
            (getline(myfile, line))
                if(line != "N")
                    return;
            N = atoi(line.c_str());

            //get A
            getline(myfile,line));
            A = new double[N][N];

            //read each ith row
            for(int i = 0; i < N; i++) {
                getline(myfile,line);
                string buf;
                stringstream ss(line);
                vector<string> row;
                
                while(ss >> buf)
                    row.push_back(buf);
                for(int j = 0; j < 
                A[i

                
            }

            myfile.close();
        }

        for (int q = 1; q < noOfProcesses; q++) {
            //MPI_Send(msg, str_len, MPI_CHAR, q, 0, MPI_COMM_WORLD);
        }

    } else {

        //MPI_Recv(msg, str_len, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    e_time = MPI_Wtime();

    if (myId != 0)
        printf("Test: process %d received message %ld bytes in %f seconds\n", myId, str_len, e_time - s_time);

    MPI_Finalize();
    return 0;
}

