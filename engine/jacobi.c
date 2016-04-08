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
	float h = 1.0 / grid->len;
	grid->value_prevs[x][y] = grid->values[x][y];
	grid->values[x][y] = 0;
	if(grid->neumann_presents[x][y] & (1 << DOWN)) {
		grid->values[x][y] += grid->value_prevs[x][y];
		grid->values[x][y] -= h * grid->neumanns[DOWN][x][y];
	} else {
		if(y < grid->len-1)
			grid->values[x][y] += grid->values[x][y+1];
	}

	if(grid->neumann_presents[x][y] & (1 << UP)) {
		grid->values[x][y] += grid->value_prevs[x][y];
		grid->values[x][y] -= h * grid->neumanns[UP][x][y];
	} else {
		if(y > 0)
			grid->values[x][y] += grid->value_prevs[x][y-1];
	}

	if(grid->neumann_presents[x][y] & (1 << RIGHT)) {
		grid->values[x][y] += grid->value_prevs[x][y];
		grid->values[x][y] -= h * grid->neumanns[RIGHT][x][y];
	} else {
		if(x < grid->len-1)
			grid->values[x][y] += grid->values[x+1][y];
	}

	if(grid->neumann_presents[x][y] & (1 << LEFT)) {
		grid->values[x][y] += grid->value_prevs[x][y];
		grid->values[x][y] -= h * grid->neumanns[LEFT][x][y];
	} else {
		if(x > 0)
			grid->values[x][y] += grid->value_prevs[x-1][y];
	}

	grid->values[x][y] += pow(h, 2.0) * grid->initials[x][y];
	grid->values[x][y] /= 4.f;
}

static inline void cell_normal(struct grid *grid, int x, int y)
{
	grid->value_prevs[x][y] = grid->values[x][y];
	grid->values[x][y] = 0.f;
	if(y < grid->len-1)
		grid->values[x][y] += grid->values[x][y+1];
	if(y > 0)
		grid->values[x][y] += grid->value_prevs[x][y-1];
	if(x < grid->len-1)
		grid->values[x][y] += grid->values[x+1][y];
	if(x > 0)
		grid->values[x][y] += grid->value_prevs[x-1][y];

	grid->values[x][y] += pow((1.0 / grid->len), 2.0) * grid->initials[x][y];
	grid->values[x][y] /= 4.f;
}

static inline double do_cell(struct grid *grid, int x, int y)
{
	if(grid->dirichlet_presents[x][y]) {
		grid->value_prevs[x][y] = grid->values[x][y] = grid->dirichlets[x][y];
		return 0.f;
	}
	if(x == 0 || y == 0 || x == grid->len-1 || y == grid->len-1)
		return 0.f;
	double old = grid->values[x][y];
	if(!grid->neumann_presents[x][y]) {
		/* do the thing! */
		cell_normal(grid, x, y);
		return pow(old - grid->values[x][y], 2.0);
	} else {
		cell_neumann(grid, x, y);
		return pow(old - grid->values[x][y], 2.0);
	}
}

static _Atomic double thresh = 0.000001;
#define THREADS 0

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

#define NUM_THREADS 4

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
			if(cur_iter++ % 100 == 0 && cur_iter > 0) {
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
			grid->values[x][y] = 0.f;
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
			grid->values[x][y] = 0.f;
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
		if((++iter % 10) == 0) {
			printf("%d: err: %.15lf\n", iter, iter_err);
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

