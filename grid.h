#pragma once

#include <stddef.h>
#include <stdbool.h>

struct cell {
	float value;
	float neumann[4];
	float value_prev;
	float dirichlet;
	float initial;
	float error;
	char neumann_present[4];
	char dirichlet_present;
};

struct grid {
	int len;
	struct cell *cells[];
};

