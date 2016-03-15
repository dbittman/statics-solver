#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "grid.h"
#define UP 0
#define RIGHT 1
#define DOWN 2
#define LEFT 3
#include <pmmintrin.h>
#include <immintrin.h>
/* do we need to store error per cell? */

struct params {
	const float omega;
	const float h;
} *params;

#include <assert.h>
static inline void cell_neumann(struct grid *grid, int x, int y)
{
#if 0
	register float value_prev = (*grid->values)[x][y];
	(*grid->values)[x][y] = 0;
	if(grid->neumann_presents) {
		(*grid->values)[x][y] += value_prev;
		(*grid->values)[x][y] -= params->h * cell->neumann[DOWN];
	} else {
		(*grid->values)[x][y] += (*grid->values)[x][y+1];
	}

	if(cell->neumann_present[UP]) {
		(*grid->values)[x][y] += value_prev;
		(*grid->values)[x][y] -= params->h * cell->neumann[UP];
	} else {
		(*grid->values)[x][y] += (*grid->values)[x][y-1];
	}

	if(cell->neumann_present[RIGHT]) {
		(*grid->values)[x][y] += value_prev;
		(*grid->values)[x][y] -= params->h * cell->neumann[RIGHT];
	} else {
		(*grid->values)[x][y] += (*grid->values)[x+1][y];
	}

	if(cell->neumann_present[LEFT]) {
		(*grid->values)[x][y] += value_prev;
		(*grid->values)[x][y] -= params->h * cell->neumann[LEFT];
	} else {
		(*grid->values)[x][y] += (*grid->values)[x-1][y];
	}

	(*grid->values)[x][y] += pow(params->h, 2.0) * cell->initial;
	(*grid->values)[x][y] *= params->omega / 4.f;
	(*grid->values)[x][y] += (1.f - params->omega) * value_prev;
#endif
}

static inline float sum(__m256 x) {
	__m128 hi = _mm256_extractf128_ps(x, 1);
	__m128 lo = _mm256_extractf128_ps(x, 0);
	lo = _mm_add_ps(hi, lo);
	hi = _mm_movehl_ps(hi, lo);
	lo = _mm_add_ps(hi, lo);
	hi = _mm_shuffle_ps(lo, lo, 1);
	lo = _mm_add_ss(hi, lo);
	return _mm_cvtss_f32(lo);
}

static inline void cell_normal(struct grid *grid, int x, int y)
{
	__mmask8 mask = 0xF;

	__attribute__((aligned(16))) float vec[8] = {
		grid->values[x][y+1],
		grid->values[x][y-1],
		grid->values[x+1][y],
		grid->values[x-1][y],
		params->h * grid->initials[x][y],
		((1.f - params->omega) / (params->omega / 4.)) * grid->values[x][y], 0., 0.
	};
	__m256 v1 = _mm256_load_ps(vec);
	__m256 scal = _mm256_set1_ps(params->omega / 4.f);
	__m256 res = _mm256_mul_ps(v1, scal);
	grid->values[x][y] = sum(res);
}

static inline double do_cell(struct grid *grid, int x, int y)
{
	if(grid->dirichlet_presents[x][y]) {
		grid->values[x][y] = grid->dirichlets[x][y];
		return 0.f;
	}
	if(x == 0 || y == 0 || x == grid->len-1 || y == grid->len-1)
		return 0.f;
	float old = grid->values[x][y];
	if(!grid->neumann_presents[x][y]) {
		/* do the thing! */
		cell_normal(grid, x, y);
		return pow(old - grid->values[x][y], 2.0);
	} else {
		cell_neumann(grid, x, y);
		return pow(old - grid->values[x][y], 2.0);
	}
}

double sor_solve(struct grid *grid)
{
	int iter = 0;
	double iter_err;
	double thresh = 0.0000001;
	for(int x=0;x<grid->len;x++) {
		for(int y=0;y<grid->len;y++) {
			grid->values[x][y] = grid->value_prevs[x][y] = 0.f;
		}
	}
	params = malloc(sizeof(struct params));
	struct params init = { .omega = 2.f / (1.f + M_PI / grid->len), .h = pow((1.f / grid->len), 2.0) };
	memcpy(params, &init, sizeof(init));
	printf("omega = %f\n", params->omega);
	do {
		iter_err = 0.f;
		for(int x=0;x<grid->len;x++) {
			for(int y=0;y<grid->len;y++) {
				iter_err += do_cell(grid, x, y);
			}
		}
		iter_err = sqrt(iter_err);
		if((++iter % 100) == 0) {
			printf("%d: err_rem: %.12lf\r", iter, iter_err - thresh);
			fflush(stdout);
			if(stop)
				break;
		}
	} while(iter_err >= thresh);
	grid->iters = iter;
	printf("\n");
	return iter_err;
}

