// mpic++ tour10.c -fopenmp -o t10
// for ip in `cat m.txt`; do scp t10 $ip:;done
// mpirun --machinefile m.txt -pernode t10 <N>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <omp.h>
#include <unistd.h>
#include <mpi.h>
using namespace std;

#define MASTER 0

struct Tabuleiro{
	int N;
  int **tab;
};

// Below arrays details all 8 possible movements for a knight.
// It is important to avoid changing sequence of below arrays
int xMove[] = {2, 1, -1, -2, -2, -1, 1, 2, 2};
int yMove[] = {1, 2, 2, 1, -1, -2, -2, -1, 1};


// Check if (x, y) is valid chess board coordinates
// Note that a knight cannot go out of the chessboard
struct Task
{
    int x;
    int y;
    int last_x;
    int last_y;
    int atual_x;
    int atual_y;
};

bool isValid(int x, int y, int N);
void printMatriz(int **visited, int N);
void printTask(Task task);

void SaidaCerta(int **sol, int n);
void SaidaCertaRot4(int **sol, int N);
void SaidaCertaRot8(int **sol, int N);
int DelegaThreads(int x, int y, int modo, int N);
Task CreateTask(int x, int y, int last_x, int last_y, int newX, int newY);
int TourSubArvore(int **visited, int N, int x, int y, int pos, int total_solucoes, int modo);
int MetaTourSubArvore(Task task, int modo, int N);

bool isValid(int x, int y, int N)
{
    if (x < 0 || y < 0 || x >= N || y >= N)
        return false;

    return true;
}

void printMatriz(int **visited, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            printf("%2d ", visited[i][j]);
        printf("\n");
    }
    printf("\n");
}

void printTask(Task task)
{
    printf("x: %d\n", task.x);
    printf("y: %d\n", task.y);
    printf("last_x: %d\n", task.last_x);
    printf("last_y: %d\n", task.last_y);
    printf("atual_x: %d\n", task.atual_x);
    printf("atual_y: %d\n\n", task.atual_y);
}

int DelegaThreads(int x, int y, int modo, int N)
{

    int **visited  = (int **)malloc(N * sizeof(int*));
    for(int i = 0; i < N; i++){
        visited[i] = (int *)malloc(N * sizeof(int));
        for(int b = 0; b < N; b++){
            visited[i][b]=-1;
        }
    }

    int total_solucoes = 0;
    vector<Task> tasks;
    visited[x][y] = 0;

    for (int k = 0; k < 8; k++)
    {
        // Get the new position of Knight from current
        // position on chessboard
        int newX1 = x + xMove[k];
        int newY1 = y + yMove[k];

        // if new position is a valid and not visited yet
        if (isValid(newX1, newY1, N) && (visited[newX1][newY1] == -1))
        {
            visited[newX1][newY1] = 1;
            for (int t = 0; t < 8; t++)
            {
                // Get the new position of tnight from current
                // position on chessboard
                int newX2 = newX1 + xMove[t];
                int newY2 = newY1 + yMove[t];

                // if new position is a valid and not visited yet
                if (isValid(newX2, newY2, N) && (visited[newX2][newY2] == -1))
                {
                    tasks.push_back(CreateTask(x, y, newX1, newY1, newX2, newY2));
                }
                visited[newX1][newY1] = -1;
            }
        }
    }

#pragma omp parallel for schedule(dynamic) reduction(+ \
                                                     : total_solucoes) num_threads(8)
    for (int i = 0; i < int(tasks.size()); i++)
    {
        total_solucoes += MetaTourSubArvore(tasks[i], modo, N);
    }
    return total_solucoes;
}

Task CreateTask(int x, int y, int last_x, int last_y, int newX, int newY)
{
    Task task;
    task.x = x;
    task.y = y;
    task.last_x = last_x;
    task.last_y = last_y;
    task.atual_x = newX;
    task.atual_y = newY;
    return task;
}

int MetaTourSubArvore(Task task, int modo, int N)
{
    int **visited  = (int **)malloc(N * sizeof(int*));
    for(int i = 0; i < N; i++){
        visited[i] = (int *)malloc(N * sizeof(int));
        for(int b = 0; b < N; b++){
            visited[i][b]=-1;
        }
    }

    visited[task.x][task.y] = 0;
    visited[task.last_x][task.last_y] = 1;
    //printMatriz(visited);
    //printf("%d %d %d %d %d %d\n",task.x, task.y, task.last_x, task.last_y, task.atual_x, task.atual_y);
    int total = TourSubArvore(visited, N, task.atual_x, task.atual_y, 2, 0, modo);
    //printf("total: %d\n", total);
    return total;
}

// Recursive function to perform the Knight's tour using backtracking
int TourSubArvore(int **visited, int N, int x, int y, int pos, int total_solucoes, int modo)
{
    // mark current square as visited

    visited[x][y] = pos;

    // if all squares are visited, print the solution
    if (pos >= (N * N) - 1)
    {
        //printMatriz(visited);
        if (modo == 4)
        {
            total_solucoes+=4;
            SaidaCerta(visited, N);
            SaidaCertaRot4(visited, N);
        }
        else if (modo == 8)
        {
            total_solucoes+=8;
            SaidaCerta(visited, N);
            SaidaCertaRot8(visited, N);
        }
        else
        {
            total_solucoes++;
            SaidaCerta(visited, N);
        }

        // backtrack before returning
        visited[x][y] = -1;
    }

    // check for all 8 possible movements for a knight
    // and recur for each valid movement
    for (int k = 0; k < 8; k++)
    {
        // Get the new position of Knight from current
        // position on chessboard
        int newX = x + xMove[k];
        int newY = y + yMove[k];

        // if new position is a valid and not visited yet
        if (isValid(newX, newY, N) && (visited[newX][newY] == -1))
            total_solucoes = TourSubArvore(visited, N, newX, newY, pos + 1, total_solucoes, modo);
    }

    // backtrack from current square and remove it from current path
    visited[x][y] = -1;
    return total_solucoes;
}

void SaidaCerta(int **sol, int N)
{
    int ajuda[N * N], k;

    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < N; y++)
        {
            ajuda[sol[x][y]] = x*N + y;
        }
    }
    #pragma omp critical
    {
         MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

void SaidaCertaRot4(int **sol, int N)
{
    int x, y, k;
    int **temp = (int **)malloc(N * sizeof(int*));
    for(int i = 0; i < N; i++)
        temp[i] = (int *)malloc(N * sizeof(int));

    int **s = (int **)malloc(N * sizeof(int*));
    for(int i = 0; i < N; i++)
        s[i] = (int *)malloc(N * sizeof(int));

    int ajuda[N*N];
    //preenche s
    for(int a = 0; a < N; a++){
        for(int b = 0; b < N; b++){
            s[a][b] = sol[a][b];
        }
    }

    //faz 3 rotacoes
    for (int i = 0; i < 3; i++)
    {
        for (x = 0; x < N; x++)
        {
            for (y = N - 1; y >= 0; y--)
            {
                temp[x][N - y - 1] = s[y][x];
                ajuda[s[y][x]] = x * N + y;
            }
        }
        #pragma omp critical
        {
             MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        for(int a = 0; a < N; a++){
            for(int b = 0; b < N; b++){
                s[a][b] = temp[a][b];
            }
        }
    }
}

void SaidaCertaRot8(int **sol, int N)
{
    int x, y, k;

    int **temp = (int **)malloc(N * sizeof(int*));
    for(int i = 0; i < N; i++)
        temp[i] = (int *)malloc(N * sizeof(int));

    int **s = (int **)malloc(N * sizeof(int*));
    for(int i = 0; i < N; i++)
        s[i] = (int *)malloc(N * sizeof(int));

    int **temp2 = (int **)malloc(N * sizeof(int*));
    for(int i = 0; i < N; i++)
        temp2[i] = (int *)malloc(N * sizeof(int));

    int **s2 = (int **)malloc(N * sizeof(int*));
    for(int i = 0; i < N; i++)
        s2[i] = (int *)malloc(N * sizeof(int));

    int ajuda[N * N];
    int ajuda2[N * N];
    //preenche s
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            s[i][j] = sol[i][j];
            s2[j][i] = sol[i][j];
        }
    }

    SaidaCerta(s2, N); // imprime inverso

    //faz 3 rotacoes da pos normal e inversa
    for (int i = 0; i < 3; i++)
    {
        for (x = 0; x < N; x++)
        {
            for (y = N - 1; y >= 0; y--)
            {
                temp[x][N - y - 1] = s[y][x];
                temp2[x][N - y - 1] = s2[y][x];
                ajuda[s[y][x]] = x * N + y;
                ajuda2[s2[y][x]] = x * N + y;
            }
        }

        #pragma omp critical
        {
             MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
             MPI_Send(ajuda2, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        for(int a = 0; a < N; a++){
            for(int b = 0; b < N; b++){
                s[a][b] = temp[a][b];
                s2[a][b] = temp2[a][b];
            }
        }

    }
}

int main(int argc, char **argv)
{
    int taskid, ntasks;
    double start_time, end_time;
    FILE *fptr;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    if ((argc != 2)) {
        printf("Uso: mpirun --machinefile m.txt -pernode executavel <tamanho do tabuleiro>\n");
        //MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }

    int N = atoi(argv[1]);

    if (N<5) {
        printf("Não há solucoes para N=%s\n", argv[0]);
        //MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }

    if (taskid == MASTER) {
        fptr = fopen("Gabriel_Marcelo-ver.txt", "w");
        fclose(fptr);
        fptr = fopen("Gabriel_Marcelo-cont.txt", "w");
        fclose(fptr);
    }


    int ajuda[N*N] = {};
    int total_solucoes = 0;
    start_time = MPI_Wtime();
    if (N == 5) {
        switch (taskid) {
            case 0:
                printf("Começou...\n");
                break;
            case 1:
                total_solucoes += DelegaThreads(1, 0, 8, N);
                total_solucoes += DelegaThreads(2, 1, 4, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 2:
                total_solucoes += DelegaThreads(1, 1, 4, N);
                total_solucoes += DelegaThreads(2, 2, 1, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 3:
                total_solucoes += DelegaThreads(2, 0, 4, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 4:
                total_solucoes += DelegaThreads(0, 0, 4, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
        }
    }else if(N == 6){
        switch (taskid) {
            case 0:
            printf("Começou...\n");
                break;
            case 1:
               total_solucoes += DelegaThreads(0, 0, 4, N);
               total_solucoes += DelegaThreads(1, 0, 8, N);
               MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 2:
               total_solucoes += DelegaThreads(1, 1, 4, N);
               MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 3:
                total_solucoes += DelegaThreads(2, 2, 4, N);
                total_solucoes += DelegaThreads(2, 1, 8, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 4:
                total_solucoes += DelegaThreads(2, 0, 8, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
        }
    }else if(N == 7){
        switch (taskid) {
            case 0:
            printf("Começou...\n");
                break;
            case 1:
                total_solucoes += DelegaThreads(0, 0, 4, N);
                total_solucoes += DelegaThreads(1, 0, 8, N);
                total_solucoes += DelegaThreads(1, 1, 4, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 2:
                total_solucoes += DelegaThreads(2, 1, 8, N);
                total_solucoes += DelegaThreads(2, 2, 4, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 3:
                total_solucoes += DelegaThreads(2, 0, 8, N);
                total_solucoes += DelegaThreads(3, 0, 4, N);
                total_solucoes += DelegaThreads(3, 1, 4, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 4:
                total_solucoes += DelegaThreads(3, 2, 4, N);
                total_solucoes += DelegaThreads(3, 3, 1, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
        }
    }else if(N == 8){
        switch (taskid) {
            case 0:
            printf("Começou...\n");
                break;
            case 1:
                total_solucoes += DelegaThreads(0, 0, 4, N);
                total_solucoes += DelegaThreads(1, 0, 8, N);
                total_solucoes += DelegaThreads(1, 1, 4, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 2:
                total_solucoes += DelegaThreads(2, 0, 8, N);
                total_solucoes += DelegaThreads(2, 1, 8, N);
                total_solucoes += DelegaThreads(2, 2, 4, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 3:
                total_solucoes += DelegaThreads(3, 0, 8, N);
                total_solucoes += DelegaThreads(3, 1, 8, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
            case 4:
                total_solucoes += DelegaThreads(3, 2, 8, N);
                total_solucoes += DelegaThreads(3, 3, 4, N);
                MPI_Send(ajuda, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
                break;
        }
    }else {
        if (taskid == MASTER) {
            printf("O algoritmo não esta configurado para calcular para N=%d\n", N);
            MPI_Finalize();
            return 0;
        }
    }

    int solsMPI = 0;
    MPI_Status status;

    int print[N*N];
    MPI_Request request2;
    if (taskid == MASTER) {
       fptr = fopen("Gabriel_Marcelo-ver.txt", "w");
       int l = 0, proc = 0;
       while (true) {
          MPI_Recv(print, N*N, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
          l++;
          if (print[0]==0 && print[1]==0) {
             proc++;
          }else{
            printf("%10d: ", l);
            for (int k = 0; k < N * N; k++)
            {
                fprintf(fptr, "%2d ", print[k]);
                printf("%2d ", print[k]);
            }
            fprintf(fptr, "\n");
            printf("\n");
          }
          if (proc>=4) {
            break;
          }
       }
       fclose(fptr);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (taskid != MASTER) {
       MPI_Send(&total_solucoes, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
       //printf("SOLUCOES %d of %d: %d\n", taskid, ntasks, solsMPI);
    } else{
        //printf("SOLUCOES %d: %d\n", taskid, total_solucoes);
        for (int i = 1; i < ntasks; i++) {
            MPI_Recv(&solsMPI, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            printf("SOLUCOES %d: %d\n", i, solsMPI);
            total_solucoes += solsMPI;
       }
    }

    end_time = MPI_Wtime();

    if (taskid == MASTER) {
       printf("\nTOTAL DE SOLUCOES %d (%d X %d): %d\nTempo: %f\n", taskid, N, N, total_solucoes, (end_time-start_time));
       fptr = fopen("Gabriel_Marcelo-cont.txt", "a");
       fprintf(fptr, "\nTOTAL DE SOLUCOES(%d X %d): %d\nTempo: %f\n", N, N, total_solucoes, (end_time-start_time));
       fclose(fptr);
    }
    MPI_Finalize();
    return 0;
}
