#include <stdio.h>
#include <stdlib.h>
#include "tsp.h"
#include "parser.h"

#define MAX_DIM 200
#define MAX_POP 256
#define TERM_NUM 200
#define BIG_NUM 10000

typedef struct {
	int c[MAX_DIM];
	double fit;
} Indiv;

typedef struct {
	int amount;
	Indiv v[2 * MAX_POP];
} Popu;

void randperm(int c[], int a, int b)
{
	int mark[MAX_DIM];
	int i, j, cnt;

	if (c == NULL || a > b) {
		printf("randperm error\n");
		exit(1);
	}

	for (i = a; i <= b; i++) {
		mark[i] = 0;
	}

	for (i = 0; i <= b - a; i++) {
		r = rand() % (b - a + 1 - i);
		cnt = 0;
		for (j = a; j <= b && cnt < r - 1; j++) {
			if (mark[j] != 1) {
				cnt++;
			}
		}
		while (mark[j] == 1) {j++;}
		c[i] = j;
	}
}

void fitSort(Indiv v[], int low, int high)	
{
	int i, j;
	Indiv tmp;

	if (low < high) {
		i = low;
		j = high;
		memcpy(&tmp, &v[i], sizeof(Indiv));
		while (i < j) {
			while (i < j && v[j].fit < tmp.fit) {
				j--;
			}
			if (i < j) {
				memcpy(&v[i], &v[j], sizeof(Indiv));
			}
			i++;
			while (i < j && v[i].fit > tmp.fit) {
				i++;
			}
			if (i < j) {
				memcpy(&v[j], &v[i], sizeof(Indiv));
			}
			j--;
		}
		memcpy(&v[i], &tmp, sizeof(Indiv));
		fitSort(v, low, i - 1);
		fitSort(v, i + 1, high);
	}
}

double getFitness(Mtx *m, int c[])
{
	int i, j;
	double fit;

	fit = 0;
	for (i = 0; i < m->dim - 1; i++) {
		fit += m->m[c[i]][c[i + 1]];
	}
	fit += m->m[c[m->dim - 1]][c[0]];

	return (1.0 / fit);
}

void initPopu(const Mtx *m, Popu *p, int popSize)
{
	int i;

	for (i = 0; i < popSize; i++) {
		randperm(p->v[i].c, 1, m->dim);
		p->v[i].fit = getFitness(m, p->v[i].c);
	}

	p->amount = popSize;
}

void speSelect(const Mtx *m, Popu *p, int maxPopSize)
{
	int i;

	//calculate children's fitness
	for (i = maxPopSize; i < p->amount; i++) {
		p->v[i].fit = getFitness(m, p->v[i].c);
	}
	
	//sort fitness
	fitSort(p->v, 0, p->amount - 1);	

	//select
	p->amount = maxPopSize;
}

void crossOver(const Mtx *m, Popu *p, int maxPopSize, double crossProb, double matingProb)
{
	
}

void speMutat(const Mtx *m, Popu *p, int maxPopSize, double mutatProb)
{
	int i, j;
	int r1, r2, t;

	for (i = maxPopSize; i < p->amount; i++) {
		if ((rand() % BIG_NUM) / BIG_NUM < mutatpRrob) {
			r1 = rand() % m->dim;
			r2 = rand() % m->dim;
			if (r1 != r2) {
				t = p->v[i].c[r1];
				p->v[i].c[r1] = p->v[i].c[r2];
				p->v[i].c[r2] = t;
			}
		}
	}
}

void showCurInfo(const Mtx *m, Popu *p, int gen)
{
	int i;

	print("Generation %d:\n", gen);
	printf("max fitness: %.3lf, min distance: %.3lf\n", p->v[0].fit, 1.0 / p->v[0].fit);
	printf("tour:\n");
	for (i = 0; i < m->dim; i++) {
		printf("%4d", p->v[0].c[i]);
		if ((i + 1) % 20 == 0) {
			printf("\n");
		}
	}
}

int main()
{ 
	Mtx mtx;
	Popu p;
	int i;

	allocMtx("eil51.tsp", &mtx);
	if (mtx->m == NULL) {
		printf("alloc matrix error\n");
		return 1;
	}

	srand(time(NULL));

	initPopu(&mtx, &p, MAX_POP);
	for (i = 0; i < TERM_NUM; i++) {
		speSelect(&mtx, &p, MAX_POP);
		showCurInfo(&mtx, &p, i);
		crossOver(&mtx, &p, MAX_POP, 0.8, 0.8);
		speMutat(&mtx, &p, MAX_POP, 0.1);
	}

	freeMtx(&mtx);
	return 0;
}
