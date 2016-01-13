#include <stdio.h>
#include <stdlib.h>
#include <stdalign.h>
#include <string.h>
#include <math.h>
int GRID_LENGTH;

#define MAX_GRID_LENGTH 5000

alignas (16) double grid[MAX_GRID_LENGTH][MAX_GRID_LENGTH][2];
alignas (16) double init[MAX_GRID_LENGTH][MAX_GRID_LENGTH];

double h;

void writeout(int iter)
{
	char name[64];
	if(iter == -1)
		sprintf(name, "out");
	else
		sprintf(name, "out%d", iter);
	FILE *f = fopen(name, "w");
	
	for(int x=0;x<GRID_LENGTH;x++) {
		for(int y=0;y<GRID_LENGTH;y++) {
			fprintf(f, "%d\t%d\t%lf\n", x, y, grid[y][x][0]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

int main(int argc, char **argv)
{

	GRID_LENGTH = atoi(argv[1]);
	h = 1.0 / (GRID_LENGTH + 1.0);

	for(int x=0;x<GRID_LENGTH;x++) {
		for(int y=0;y<GRID_LENGTH;y++) {
			grid[x][y][0] = grid[x][y][1] = 0.0;
			init[x][y]=0.0;
		}
	}
	for(int x=0;x<GRID_LENGTH;x++) {
		init[0][x] = 10 * sin(M_PI * (x / (double)GRID_LENGTH));
	}

	init[300][300]=40000;
	init[200][200]=40000;
	init[250][250]=-40000;
	for(int x=0;x<GRID_LENGTH;x++) {
		for(int y=0;y<GRID_LENGTH;y++) {
			grid[x][y][0] = grid[x][y][1] = init[x][y];
		}
	}


	printf("h = %lf\n", h);
	int file = 1;
	int num = 1000;
	for(int i=0;i<=num;i++) {
		fprintf(stderr, "iteration %d\r", i);
		for(int x=1;x<GRID_LENGTH-1;x++) {
			for(int y=1;y<GRID_LENGTH-1;y++) {
				grid[y][x][1] = grid[y][x][0];
				grid[y][x][0] = 0;
				int n=0;
				if(y < GRID_LENGTH-1 && x < GRID_LENGTH-1)
					grid[y][x][0] += grid[y+1][x+1][1], n++;
				if(y < GRID_LENGTH-1 && x > 0)
					grid[y][x][0] += grid[y+1][x-1][1], n++;
				if(y > 0 && x < GRID_LENGTH-1)
					grid[y][x][0] += grid[y-1][x+1][1], n++;
				if(y > 0 && x > 0)
					grid[y][x][0] += grid[y-1][x-1][1], n++;
				if(y < GRID_LENGTH-1)
					grid[y][x][0] += grid[y+1][x][1], n++;
				if(y > 0)
					grid[y][x][0] += grid[y-1][x][1], n++;
				if(x < GRID_LENGTH-1)
					grid[y][x][0] += grid[y][x+1][1], n++;
				if(x > 0)
					grid[y][x][0] += grid[y][x-1][1], n++;

				grid[y][x][0] += pow(h, 2.0) * init[y][x];
				grid[y][x][0] /= (float)n;
			}
		}

		//if(i % (num / 100) == 0)
		//	writeout(file++);
	}


	printf("\nDone\n");
	
	writeout(-1);

}

