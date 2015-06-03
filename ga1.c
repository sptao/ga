#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define POP_SIZE 128
#define CROSS_PROB 0.6
#define MATING_PROB 0.8
#define MUTAT_PROB 0.1
#define TERM_GEN 10240
#define PI 3.1415
#define DELTA 0.1
#define NUM 10000

typedef struct {
	double val;
	int idx;
} Fit;

int amount;
int child;
double chromo[POP_SIZE * 4];
Fit fit[POP_SIZE * 4];

void qsortFit(int low, int high)
{
	int i, j;
	Fit tmp;
	
	if (low < high) {
		i = low + 1;
		j = high;
		while (i < j) {
			while (j > low && fit[j].val < fit[low].val) {j--;}
			while (i < high && fit[i].val > fit[low].val) {i++;}
			if (i != j) {
				tmp = fit[j];
				fit[j] = fit[i];
				fit[i] = tmp;
			}
		}
		tmp = fit[i];
		fit[i] = fit[low];
		fit[low] = tmp;
		qsortFit(low, i - 1);
		qsortFit(i + 1, high);
	}
}

void initPop()
{
	int i;

#pragma omp unroll parallel for
	for (i = 0; i < POP_SIZE; i++) {
		chromo[i] = (rand() % NUM) / (double)NUM * 3.0 - 1.0;	//[-1, 2]
	}
	amount = POP_SIZE;
	//calculate first generation's fitness
#pragma omp unroll parallel for
	for (i = 0; i < amount; i++) {
		fit[i].val = chromo[i] * sin(10.0 * PI * chromo[i]) + 2.0;
		fit[i].idx = i;
	}
}

void calcFitness(int n)
{
	int i, j, m;
	double tc[POP_SIZE * 4];

	//only calculate child's fitness
#pragma omp unroll parallel for
	for (i = POP_SIZE; i < amount; i++) {
		fit[i].val = chromo[i] * sin(10.0 * PI * chromo[i]) + 2.0;
		fit[i].idx = i;
	}

	//sort fit according to fit[i].v
/*
	Fit tmp;
	for (i = 0; i < amount; i++) {
		m = i;
		for (j = i + 1; j < amount; j++) {
			if (fit[j].val > fit[m].val) {
				m = j;
			}
		}
		if (m != i) {
			tmp = fit[i];
			fit[i] = fit[m];
			fit[m] = tmp;
		}
	}
*/
	qsortFit(0, amount - 1);

	//select chromosome according to fit[i].v
	for (i = 0; i < POP_SIZE; i++) {
		tc[i] = chromo[fit[i].idx];
		fit[i].idx = i;
	}
	for (i = 0; i < POP_SIZE; i++) {
		chromo[i] = tc[i];
	}
	amount = POP_SIZE;

	//printf("best result after %d generation:\nx = %.6lf, f = %.6lf\n", n, chromo[0], fit[0].val);
}

void crossOver()
{
	int i, j, r;
	static int s[NUM * 2];
	int a, b, total;
	int mom[POP_SIZE];
	int dad[POP_SIZE];

	total = 0;
	for (i = 0; i < POP_SIZE; i++) {
		total += fit[i].val;
	}

	//calculate select prob
	a = 0;
	for (i = 0; i < POP_SIZE; i++) {
		b = fit[i].val / total * (double)NUM + a;
		for (j = a; j < b; j++) {
			s[j] = i;
		}
		a = b;
	}

	//select mom and dad
#pragma omp unroll parallel for
	for (i = 0; i < POP_SIZE * CROSS_PROB; i++) {
		//select mom
		mom[i] = s[rand() % NUM];
		//select dad
		dad[i] = s[rand() % NUM];
	}

#pragma omp unroll parallel for
	for (i = 0; i < POP_SIZE * CROSS_PROB; i++) {
		if ((rand() % NUM) / (double)NUM < MATING_PROB) {
			chromo[amount++] = (chromo[mom[i]] + chromo[dad[i]]) / 2.0;
		}
	}
}

void mutation(int n)
{
	int i;

	//mutation probability = 0.1
	for (i = POP_SIZE; i < amount; i++) {
		if ((rand() % NUM) / (double)NUM < MUTAT_PROB) {
			if (rand() % 2 == 0) {
				chromo[i] += DELTA / (double)(n + 1);
			}	
			else {
				chromo[i] -= DELTA / (double)(n + 1);
			}
		}
	}
}

int main()
{
	int i;

	srand(time(NULL));
	initPop();

	for (i = 0; i < TERM_GEN; i++) {
		//printf("start generation %d\n", i);
		calcFitness(i); 
		crossOver();
		mutation(i);
	}

	printf("best result: x = %.6lf, f = %.6lf\n", chromo[0], fit[0].val);

	return 0;
}
