#pragma once

#define SOLVE_METHOD_JACOBI 0
#define SOLVE_METHOD_SOR 1

double solve(struct grid *, int);
void init_grid(struct grid *);

