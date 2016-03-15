#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include "grid.h"
#define UP 0
#define RIGHT 1
#define DOWN 2
#define LEFT 3

#include <assert.h>
/* do we need to store error per cell? */
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
			cell->value += grid->cells[x][y-1].value_prev;
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
			cell->value += grid->cells[x-1][y].value_prev;
	}

	cell->value += pow(h, 2.0) * cell->initial;
	cell->value /= 4.f;
}

static inline void cell_normal(struct grid *grid, int x, int y)
{
	struct cell *cell = &grid->cells[x][y];
	cell->value_prev = cell->value;
	cell->value = 0.f;
	if(y < grid->len-1)
		cell->value += grid->cells[x][y+1].value;
	if(y > 0)
		cell->value += grid->cells[x][y-1].value_prev;
	if(x < grid->len-1)
		cell->value += grid->cells[x+1][y].value;
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
		return 0.f;
	}
	if(x == 0 || y == 0 || x == grid->len-1 || y == grid->len-1)
		return 0.f;
	bool any_neumann = false;
	for(int i=0;i<4;++i) {
		any_neumann = any_neumann || cell->neumann_present[i];
	}
	if(!any_neumann) {
		/* do the thing! */
		cell_normal(grid, x, y);
		return pow(old - grid->cells[x][y].value, 2.0);
	} else {
		cell_neumann(grid, x, y);
		return pow(old - grid->cells[x][y].value, 2.0);
	}
}

_Atomic double thresh = 0.00001;
#define THREADS 12

#if THREADS

#define _GNU_SOURCE
#include <pthread.h>

struct args {
	struct grid *grid;
	int pos;
	int sx, ex, sy, ey;
};
#include <stdatomic.h>
double cur_err, iter_err;
_Atomic int cur_iter;
_Atomic int ptit;

#define NUM_THREADS 1

void *thread_main(void *arg)
{
	struct args *args = arg;
	while(stop != true) {
		//fprintf(stderr, "Doing pos %d (%d cells)\n", args->pos, (args->ex - args->sx) * (args->ey - args->sy));
		for(int x=args->sx;x<args->ex;x++) {
			for(int y=args->sy;y<args->ey;y++) {
				cur_err += do_cell(args->grid, x, y);
			}
		}
		if(ptit++ % NUM_THREADS == 0) {
			
			//fprintf(stderr, "pos %d reset\n", args->pos);
			iter_err = sqrt(cur_err);// / pow(args->grid->len, 2.0);
			if(cur_iter++ % 100 == 0) {
				printf("%d: err: %.15lf\r", cur_iter, iter_err);
				fflush(stdout);
				if(iter_err < thresh) {
					stop = true;
					return NULL;
				}

			}
			cur_err = 0.f;
		}

	}
	return NULL;
}

double jacobi_solve(struct grid *grid)
{
	printf("threadcount = %d\n", NUM_THREADS);
	cur_err = 0.f;
	iter_err = 0.f;
	cur_iter = 0;
	ptit=0;
	for(int x=0;x<grid->len;x++) {
		for(int y=0;y<grid->len;y++) {
			grid->cells[x][y].value = 0.f;
		}
	}

	struct args args[NUM_THREADS];

	for(int i=0;i<NUM_THREADS;i++) {
		args[i].pos = i;
		args[i].grid = grid;

		args[i].sx = 0;
		args[i].ex = grid->len;
		args[i].sy = i*grid->len / NUM_THREADS;
		args[i].ey = (i+1) * grid->len / NUM_THREADS;
	}

	pthread_t threads[NUM_THREADS];
	for(int i=0;i<NUM_THREADS;i++) {
		pthread_create(&threads[i], NULL, thread_main, &args[i]);
	}

	for(int i=0;i<NUM_THREADS;i++) {
		pthread_join(threads[i], NULL);
	}
	printf("\n");
	grid->iters = cur_iter;
	return iter_err;
}
#else
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
		//fprintf(stderr, "Doing %d cells\n", grid->len * grid->len);
		for(int x=0;x<grid->len;x++) {
			for(int y=0;y<grid->len;y++) {
				iter_err += do_cell(grid, x, y);
			}
		}
		//iter_err /= pow(grid->len, 2.0);
		if((++iter % 100) == 0) {
			printf("%d: err: %.15lf\r", iter, iter_err);
			fflush(stdout);
			if(stop)
				break;
		}
	} while(iter_err >= thresh);
	grid->iters = iter;
	printf("\n");
	return iter_err;
}
#endif

