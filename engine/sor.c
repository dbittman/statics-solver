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
#define VECTORIZE 0
#define THREADS 1
#define NUM_THREADS 4
struct params {
	const float omega, prefix1, prefix2;
	const float h, h_square;
} *params;

#include <assert.h>
static inline void cell_neumann(struct grid *grid, int x, int y)
{
	register float value = 0;
	register float value_prev = grid->values[x][y];
#if VECTORIZE
	__m128 n = _mm_set_ps(
			grid->neumanns[LEFT][x][y],
			grid->neumanns[DOWN][x][y],
			grid->neumanns[RIGHT][x][y],
			grid->neumanns[UP][x][y]
			);

	__m128 scale = _mm_set1_ps(params->h);
	__m128 add = _mm_set1_ps(value_prev);

	__m128 neumanns = _mm_sub_ps(add, _mm_mul_ps(n, scale));

	__m128 values = _mm_set_ps(
			grid->values[x-1][y],
			grid->values[x][y+1],
			grid->values[x+1][y],
			grid->values[x][y-1]);

	__m128i maski = _mm_set1_epi32(grid->neumann_presents[x][y]);
	__m128 mask = _mm_castsi128_ps(-_mm_and_si128(maski, _mm_set_epi32(1 << LEFT, 1 << DOWN, 1 << RIGHT, 1 << UP)));
	
	values = _mm_blendv_ps(values, neumanns, mask);

	__m128 additionals = _mm_set_ps(0., 0., 
			params->h_square * grid->initials[x][y],
			params->prefix1 * value_prev);

	values = _mm_add_ps(values, additionals);
	values = _mm_mul_ps(values, _mm_set1_ps(params->prefix2));
	__m128 result = _mm_hadd_ps(values, values);
	value = _mm_hadd_ps(result, result)[0];
	grid->values[x][y] = value;
#else

	if(grid->neumann_presents[x][y] & (1 << DOWN)) {
		value = value_prev - params->h * grid->neumanns[DOWN][x][y];
	} else {
		value = grid->values[x][y+1];
	}

	if(grid->neumann_presents[x][y] & (1 << UP)) {
		value += value_prev;
		value -= params->h * grid->neumanns[UP][x][y];
	} else {
		value += grid->values[x][y-1];
	}

	if(grid->neumann_presents[x][y] & (1 << RIGHT)) {
		value += value_prev;
		value -= params->h * grid->neumanns[RIGHT][x][y];
	} else {
		value += grid->values[x+1][y];
	}

	if(grid->neumann_presents[x][y] & (1 << LEFT)) {
		value += value_prev;
		value -= params->h * grid->neumanns[LEFT][x][y];
	} else {
		value += grid->values[x-1][y];
	}

	value += params->h_square * grid->initials[x][y] + params->prefix1 * value_prev;
	grid->values[x][y] = value * params->prefix2;
#endif
}
#if VECTORIZE
static inline float sum8(__m256 x) {
    // hiQuad = ( x7, x6, x5, x4 )
    const __m128 hiQuad = _mm256_extractf128_ps(x, 1);
    // loQuad = ( x3, x2, x1, x0 )
    const __m128 loQuad = _mm256_castps256_ps128(x);
    // sumQuad = ( x3 + x7, x2 + x6, x1 + x5, x0 + x4 )
    const __m128 sumQuad = _mm_add_ps(loQuad, hiQuad);
    // loDual = ( -, -, x1 + x5, x0 + x4 )
    const __m128 loDual = sumQuad;
    // hiDual = ( -, -, x3 + x7, x2 + x6 )
    const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);
    // sumDual = ( -, -, x1 + x3 + x5 + x7, x0 + x2 + x4 + x6 )
    const __m128 sumDual = _mm_add_ps(loDual, hiDual);
    // lo = ( -, -, -, x0 + x2 + x4 + x6 )
    const __m128 lo = sumDual;
    // hi = ( -, -, -, x1 + x3 + x5 + x7 )
    const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);
    // sum = ( -, -, -, x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 )
    const __m128 sum = _mm_add_ss(lo, hi);
    return _mm_cvtss_f32(sum);
}
#endif
static inline void cell_normal(struct grid *grid, int x, int y)
{
#if VECTORIZE
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
	grid->values[x][y] = sum8(res);
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

void init_params(struct grid *grid)
{
	params = malloc(sizeof(struct params));
	float omega = 2.f / (1.f + M_PI / grid->len);
	struct params init = { 
		.omega = omega,
		.h = (1.f / grid->len),
		.h_square = pow((1.f / grid->len), 2.0),
		.prefix2 = omega / 4.,
		.prefix1 = (1.f - omega) / (omega / 4.),
	};
	memcpy(params, &init, sizeof(init));
	printf("omega = %f\n", params->omega);
}


static const double thresh = 0.0000000003;

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


void *thread_main(void *arg)
{
	int thisiter = 0;
	struct args *args = arg;
	while(stop != true) {
		for(int x=args->sx;x<args->ex;x++) {
			for(int y=args->sy;y<args->ey;y++) {
				cur_err += do_cell(args->grid, x, y);
			}
		}
		thisiter++;
		cur_iter++;
		if(thisiter % 1000 == 0 && args->pos == 0) {
			iter_err = sqrt(cur_err);// / pow(args->grid->len, 2.0);
			//printf("%d: err: %.15lf\r", cur_iter, iter_err);
			//fflush(stdout);
			if(cur_iter > 100000 || iter_err < thresh) {
				stop = true;
				return NULL;
			}
			cur_err = 0.f;
		}
		if(thisiter % 100 == 0)
			if(stop == true)
				break;
	}
	return NULL;
}

double sor_solve(struct grid *grid)
{
	printf("threadcount = %d\n", NUM_THREADS);
	init_params(grid);
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
	grid->iters = cur_iter / NUM_THREADS;
	return iter_err;
}
#else
double sor_solve(struct grid *grid)
{
	init_params(grid);
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
		if((++iter % 100) == 0) {
			printf("%d: err: %.15lf\n", iter, iter_err);
			if(stop || iter_err < thresh || iter >= 10000)
				break;
		}
	} while(true);
	grid->iters = iter;
	printf("\n");
	return iter_err;
}
#endif

