#pragma once

#include <stddef.h>
#include <stdbool.h>

extern _Atomic bool stop;

struct cell {
	float value;
	float value_prev;
	float dirichlet;
	float initial;
	char dirichlet_present;
	char neumann_present[4];
	float neumann[4];
};

struct grid {
	int len;
	int iters;
	struct cell *cells[];
};

