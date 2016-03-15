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
	const float omega, prefix1, prefix2;
	const float h;
} *params;

#include <assert.h>
static inline void cell_neumann(struct grid *grid, int x, int y)
{
	register float value_prev = grid->values[x][y];
	grid->values[x][y] = 0;
	if(grid->neumann_presents[x][y] & (1 << DOWN)) {
		grid->values[x][y] += value_prev;
		grid->values[x][y] -= params->h * grid->neumanns[DOWN][x][y];
	} else {
		grid->values[x][y] += grid->values[x][y+1];
	}

	if(grid->neumann_presents[x][y] & (1 << UP)) {
		grid->values[x][y] += value_prev;
		grid->values[x][y] -= params->h * grid->neumanns[UP][x][y];
	} else {
		grid->values[x][y] += grid->values[x][y-1];
	}

	if(grid->neumann_presents[x][y] & (1 << RIGHT)) {
		grid->values[x][y] += value_prev;
		grid->values[x][y] -= params->h * grid->neumanns[RIGHT][x][y];
	} else {
		grid->values[x][y] += grid->values[x+1][y];
	}

	if(grid->neumann_presents[x][y] & (1 << LEFT)) {
		grid->values[x][y] += value_prev;
		grid->values[x][y] -= params->h * grid->neumanns[LEFT][x][y];
	} else {
		grid->values[x][y] += grid->values[x-1][y];
	}

	grid->values[x][y] += pow(params->h, 2.0) * grid->initials[x][y];
	grid->values[x][y] *= params->omega / 4.f;
	grid->values[x][y] += (1.f - params->omega) * value_prev;
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
#if 0
	__m256 v1 = _mm256_set_ps(
		grid->values[x][y+1],
		grid->values[x][y-1],
		params->h * grid->initials[x][y],
		grid->values[x-1][y],
		params->prefix1 * grid->values[x][y],
		grid->values[x+1][y],
		0., 0.);
	
	__m256 scal = _mm256_set1_ps(params->prefix2);
	__m256 res = _mm256_mul_ps(v1, scal);
	grid->values[x][y] = sum(res);
#else
	float value = grid->values[x][y+1];
	value += grid->values[x][y-1];
	value += grid->values[x-1][y];
	value += grid->values[x+1][y];
	value += params->h * grid->initials[x][y];
	value *= params->omega / 4.;
	grid->values[x][y] = value + (1.f - params->omega) * grid->values[x][y];
#endif
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
	float omega = 2.f / (1.f + M_PI / grid->len);
	struct params init = { 
		.omega = omega,
		.h = pow((1.f / grid->len), 2.0),
		.prefix2 = omega / 4.,
		.prefix1 = (1.f - omega) / (omega / 4.),
	};
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

