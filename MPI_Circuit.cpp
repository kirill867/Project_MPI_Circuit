#include <stdlib.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>
#define _REENTRANT
#define r 1000
#define C 1e-6
#define T_MAX 0.0001
#define delta_t 0.00000001
#define num_I 175
int main(int argc, char **argv)
{
	int myrank, total;// идентификатор процесса и общее количество 
	MPI_Init (&argc, &argv);
	int NumberOfNodes = atoi(argv[1]);
	// массив значений – глобальный и местный (для каждого процесса)
	double *prev_u, *u; double *prev_U, *U;
	double prev, final, delta_T, U_const;
	FILE* file;
	int start_time, end_time;
	int n, m, num_all;
	MPI_Status status;
	int i;
	int intBuf[2];
	MPI_Comm_size(MPI_COMM_WORLD, &total);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (!myrank) // для root-процесса
	{
		file = fopen("gnuplot_par.plt", "w+");
		num_all = NumberOfNodes*total;
		U = (double *)malloc(sizeof(double)*num_all);
		intBuf[0] = num_all;
		intBuf[1] = NumberOfNodes;
	}
	MPI_Bcast((void *)intBuf, 2, MPI_INT, 0, MPI_COMM_WORLD);// 
	передача информации
	n = intBuf[0];
	m = intBuf[1];
	prev_u = (double *)malloc(sizeof(double)*(m));
	u = (double *)malloc(sizeof(double)*(m));
	prev = final = 0;
	for (i = 0; i<m; i++)
	{
		prev_u[i] = 0;
	}
	int NextProcNum = (myrank == total - 1) ? MPI_PROC_NULL : myrank + 1;
	int PrevProcNum = (myrank == 0) ? MPI_PROC_NULL : myrank - 1;
	if (!myrank) start_time = clock();
	int j = 0;
	for (delta_T = delta_t; delta_T <= T_MAX; delta_T += delta_t, j++)
	{
		i = 0;
		U_const = sin(2 * 3.1415 * 100000 * delta_T);
		u[i] = prev_u[i] - ((2 * prev_u[i]) / r)*(delta_T / C);
		if (!myrank)
		u[i] = 100;
		else u[i] += ((prev + prev_u[i + 1]) / r)*(delta_T / C);
		for (i = 1; i < m - 1; i++){
			float I_const = sin(2 * 3.14 * 100000 * 1);//delta_T);
			if (myrank == total / 2 && i == m / 2)
			u[i] = ((prev_u[i - 1] - prev_u[i]) / r + I_const * 10)*(delta_t / C) + prev_u[i];
			else
			u[i] = ((prev_u[i - 1] - prev_u[i]) / r - (prev_u[i] - prev_u[i + 1]) / r)*(delta_T / C) + prev_u[i];
			}
			u[i] = prev_u[i] - ((2 * prev_u[i]) / r)*(delta_T / C);
			if (myrank == total - 1){
				float I_const = sin(2 * 3.14 * 100000 * 1);
				u[i] = ((prev_u[i - 1] - prev_u[i]) / r + I_const * 
				10)*(delta_t / C) + prev_u[i];
			}
			else u[i] += ((prev_u[i - 1] + final) / r)*(delta_T / C);
				fflush(stdout);
				if (myrank != total - 1)
				MPI_Sendrecv(&u[m - 1], 1, MPI_DOUBLE, NextProcNum, 2,
				&final, 1, MPI_DOUBLE, NextProcNum, 2, MPI_COMM_WORLD, 
				&status);
				fflush(stdout);
			if (myrank != 0)
				MPI_Sendrecv(&u[0], 1, MPI_DOUBLE, PrevProcNum, 2,
				&prev, 1, MPI_DOUBLE, PrevProcNum, 2, MPI_COMM_WORLD, 
				&status);
				MPI_Gather((void *)u, m, MPI_DOUBLE, (void *)U, m, MPI_DOUBLE, 
				0, MPI_COMM_WORLD); //сбор информации
			if (!myrank && (j % 10 == 0))
			{
			for (int i = 0; i < n; i++)
				fprintf(file, "%d %.10f\n", i, U[i]);
				fprintf(file, "\n\n");
			}
			for (i = 0; i < m; i++)
				prev_u[i] = u[i];
			}
			if (!myrank)
			{
				end_time = clock();
				printf("Processing time: %d\n", end_time - start_time);
			}
			if (!myrank)
			{
				FILE *fid = fopen("comm", "w+"); // открываем файл для записи
				формул gnuplot
				fprintf(fid, "set xrange[0:%d]\n", n - 1);
				fprintf(fid, "set yrange[-200:200]\n");
				fprintf(fid, "do for [i=0:%d]{\n", (int)(T_MAX / (10 * delta_t)));
				fprintf(fid, "plot 'gnuplot_par.plt' index i using 1:2 smooth 
				bezier\npause 0.1}\npause -1\n"); fclose(fid);
			}
			if (!myrank)
				fclose(file);
				MPI_Finalize();// завершаем работу MPI
				exit(0);
}