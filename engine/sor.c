#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include "grid.h"
#define UP 0
#define RIGHT 1
#define DOWN 2
#define LEFT 3

/* do we need to store error per cell? */

static float omega = 1.5;

static inline void cell_neumann(struct grid *grid, int x, int y)
{
	struct cell *cell = &grid->cells[x][y];
	float h = 1.0 / grid->len;
	cell->value_prev = cell->value;
	cell->value = 0;
	if(cell->neumann_present[DOWN]) {
		cell->value += cell->value_prev;
		cell->value -= h * cell->neumann[DOWN];
	} else {
		if(y < grid->len-1)
			cell->value += grid->cells[x][y+1].value;
	}

	if(cell->neumann_present[UP]) {
		cell->value += cell->value_prev;
		cell->value -= h * cell->neumann[UP];
	} else {
		if(y > 0)
			cell->value += grid->cells[x][y-1].value;
	}

	if(cell->neumann_present[RIGHT]) {
		cell->value += cell->value_prev;
		cell->value -= h * cell->neumann[RIGHT];
	} else {
		if(x < grid->len-1)
			cell->value += grid->cells[x+1][y].value;
	}

	if(cell->neumann_present[LEFT]) {
		cell->value += cell->value_prev;
		cell->value -= h * cell->neumann[LEFT];
	} else {
		if(x > 0)
			cell->value += grid->cells[x-1][y].value;
	}

	cell->value += pow(h, 2.0) * cell->initial;
	cell->value *= omega / 4.f;
	cell->value += (1.f - omega) * cell->value_prev;
}

static inline void cell_normal(struct grid *grid, int x, int y)
{
	struct cell *cell = &grid->cells[x][y];
	cell->value_prev = cell->value;
	cell->value = 0;
	if(y < grid->len-1)
		cell->value += grid->cells[x][y+1].value;
	if(y > 0)
		cell->value += grid->cells[x][y-1].value;
	if(x < grid->len-1)
		cell->value += grid->cells[x+1][y].value;
	if(x > 0)
		cell->value += grid->cells[x-1][y].value;

	cell->value += pow((1.0 / grid->len), 2.0) * cell->initial;
	cell->value *= omega / 4.f;
	cell->value += (1.f - omega) * cell->value_prev;
}

static inline double do_cell(struct grid *grid, int x, int y)
{
	struct cell *cell = &grid->cells[x][y];
	double old = grid->cells[x][y].value;
	if(cell->dirichlet_present) {
		cell->value_prev = cell->value = cell->dirichlet;
		cell->error = 0.f;
		return 0.f;
	}
	bool any_neumann = false;
	for(int i=0;i<4;++i) {
		any_neumann = any_neumann || cell->neumann_present[i];
	}
	if(!any_neumann) {
		/* do the thing! */
		cell_normal(grid, x, y);
		cell->error = fabs(old - grid->cells[x][y].value);
		return cell->error;
	} else {
		cell_neumann(grid, x, y);
		cell->error = fabs(old - grid->cells[x][y].value);
		return cell->error;
	}
}

double sor_solve(struct grid *grid)
{
	int iter = 0;
	double iter_err;
	double thresh = 0.000001;
	for(int x=0;x<grid->len;x++) {
		for(int y=0;y<grid->len;y++) {
			grid->cells[x][y].value = grid->cells[x][y].value_prev = 0.f;
		}
	}
	omega = 2.f / (1.f + M_PI / grid->len);
	printf("omega = %f\n", omega);
	do {
		iter_err = 0.f;
		for(int x=0;x<grid->len;x++) {
			for(int y=0;y<grid->len;y++) {
				iter_err += do_cell(grid, x, y);
			}
		}
		iter_err /= pow(grid->len, 2.0);
		if((++iter % 100) == 0) {
			printf("%d: err_rem: %.12lf\r", iter, iter_err - thresh);
			fflush(stdout);
		}
	} while(iter_err >= thresh);
	grid->iters = iter;
	printf("\n");
	return iter_err;
}

