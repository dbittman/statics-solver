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

