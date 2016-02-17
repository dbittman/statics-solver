#include <stdio.h>
#include <stdlib.h>
#include <stdalign.h>
#include <string.h>
#include <math.h>
#include "grid.h"

extern double jacobi_solve(struct grid *grid);
extern double sor_solve(struct grid *grid);

const char *solvers_names[] = {
	"jacobi",
	"SOR",
	"FFT",
	"magic",
};

static double (*solvers[])(struct grid *) = {
	jacobi_solve,
	sor_solve,
	NULL,
	NULL
};

double solve(struct grid *grid, int method)
{
	printf("Started SOLVER - gridsz %d, method %d\n", grid->len, method);
	if(solvers[method] == NULL) {
		fprintf(stderr, "Solver method %d (%s) not implemented\n", method, solvers_names[method]);
		return NAN;
	}
	double res = solvers[method](grid);
	printf("SOLVER completed after %d iterations with err %lf\n", grid->iters, res);
	return res;
}

