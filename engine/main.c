#include <stdio.h>
#include <stdlib.h>
#include <stdalign.h>
#include <string.h>
#include <math.h>
#include "grid.h"
#include <signal.h>
extern double jacobi_solve(struct grid *grid);
extern double sor_solve(struct grid *grid);

_Atomic bool stop = false;

const char *solvers_names[] = {
	"jacobi",
	"SOR",
	"FFT",
	"magic",
};

void sigint(int sig)
{
	stop = true;
}

static double (*solvers[])(struct grid *) = {
	jacobi_solve,
	sor_solve,
	NULL,
	NULL
};
#include <sys/time.h>
double solve(struct grid *grid, int method)
{
	struct timeval t1, t2;
    double elapsedTime;
    gettimeofday(&t1, NULL);
	signal(SIGINT, sigint);
	printf("Started SOLVER - gridsz %d, method %d\n", grid->len, method);
	if(solvers[method] == NULL) {
		fprintf(stderr, "Solver method %d (%s) not implemented\n", method, solvers_names[method]);
		return NAN;
	}
	double res = solvers[method](grid);
    gettimeofday(&t2, NULL);
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
	printf("SOLVER completed after %d iterations with err %lf, in %3.3f seconds (%f iters/sec)\n", grid->iters, res, elapsedTime / 1000.f, grid->iters / (elapsedTime / 1000.f));
	return res;
}

static void init_float(int len, float ***arr)
{
	*arr = malloc(len * sizeof(float *));
	float *tmp = malloc(len * len * sizeof(float));
	for(int i=0;i<len;i++) {
		(*arr)[i] = tmp + (i * len);
		memset((*arr)[i], 0, len * sizeof(float));
	}
}

static void init_byte(int len, uint8_t ***arr)
{
	*arr = malloc(len * sizeof(uint8_t *));
	uint8_t *tmp = malloc(len * len * sizeof(uint8_t));
	for(int i=0;i<len;i++) {
		(*arr)[i] = tmp + (i * len);
		memset((*arr)[i], 0, len * sizeof(uint8_t));
	}
}

void init_grid(struct grid *grid)
{
	printf("Initialize grid: %d\n", grid->len);
	init_float(grid->len, &grid->values);
	init_float(grid->len, &grid->value_prevs);
	init_float(grid->len, &grid->initials);
	init_float(grid->len, &grid->neumanns[0]);
	init_float(grid->len, &grid->neumanns[1]);
	init_float(grid->len, &grid->neumanns[2]);
	init_float(grid->len, &grid->neumanns[3]);
	init_float(grid->len, &grid->dirichlets);
	init_byte(grid->len, &grid->dirichlet_presents);
	init_byte(grid->len, &grid->neumann_presents);
}

