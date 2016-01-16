#pragma once

#include <stddef.h>
#include <stdbool.h>

struct cell {
	float value;
	float value_prev;
	float dirichlet;
	float initial;
	float error;
	char dirichlet_present;
	char neumann_present[4];
	float neumann[4];
};

struct grid {
	int len;
	struct cell *cells[];
};

