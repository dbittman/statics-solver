#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include "grid.h"
#define UP 0
#define RIGHT 1
#define DOWN 2
#define LEFT 3

/* do we need to store error per cell? */

static inline void cell_neumann(struct grid *grid, int x, int y)
{
	struct cell *cell = &grid->cells[x][y];
	cell->value_prev = cell->value;
	cell->value = 0.f;
	if(y < grid->len-1)
		cell->value += grid->cells[x][y+1].value_prev;
	if(y > 0)
		cell->value += grid->cells[x][y-1].value_prev;
	if(x < grid->len-1)
		cell->value += grid->cells[x+1][y].value_prev;
	if(x > 0)
		cell->value += grid->cells[x-1][y].value_prev;

	cell->value += pow((1.0 / grid->len), 2.0) * cell->initial;
	cell->value /= 4.f;
}


static inline void cell_normal(struct grid *grid, int x, int y)
{
	struct cell *cell = &grid->cells[x][y];
	cell->value_prev = cell->value;
	cell->value = 0.f;
	if(y < grid->len-1)
		cell->value += grid->cells[x][y+1].value_prev;
	if(y > 0)
		cell->value += grid->cells[x][y-1].value_prev;
	if(x < grid->len-1)
		cell->value += grid->cells[x+1][y].value_prev;
	if(x > 0)
		cell->value += grid->cells[x-1][y].value_prev;

	cell->value += pow((1.0 / grid->len), 2.0) * cell->initial;
	cell->value /= 4.f;
}

static inline double do_cell(struct grid *grid, int x, int y)
{
	struct cell *cell = &grid->cells[x][y];
	double old = grid->cells[x][y].value;
	if(cell->dirichlet_present) {
		cell->value_prev = cell->value = cell->dirichlet;
		cell->error = 0;
		return 0.f;
	}
	bool any_neumann = false;
	for(int i=0;i<4;++i) {
		any_neumann = any_neumann || cell->neumann[i];
	}
	if(!any_neumann) {
		/* do the thing! */
		cell_normal(grid, x, y);
		cell->error = fabs(old - grid->cells[x][y].value);
		return cell->error;
	}
	return 0.f;
}

double jacobi_solve(struct grid *grid)
{
	int iter = 0;
	double iter_err;
	for(int x=0;x<grid->len;x++) {
		for(int y=0;y<grid->len;y++) {
			grid->cells[x][y].value = 0.f;
		}
	}
	do {
		iter_err = 0.f;
		
		for(int x=0;x<grid->len;x++) {
			for(int y=0;y<grid->len;y++) {
				iter_err += do_cell(grid, x, y);
			}
		}
		
		iter_err /= pow(grid->len, 2.0);
		if((++iter % 100) == 0) {
			printf("%d: err: %.15lf\r", iter, iter_err);
			fflush(stdout);
		}
	} while(iter_err >= 0.0000001);
	grid->iters = iter;
	printf("\n");
	return iter_err;
}

