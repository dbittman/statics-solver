#pragma once

#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
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
	float **values;
	float **value_prevs;
	float **initials;
	uint8_t **dirichlet_presents;
	float **dirichlets;
	uint8_t **neumann_presents;
	float **neumanns[4];
};

