#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tsp.h"
#include "parser.h"

#define MAX_DIM 200
#define MAX_POP 512
#define TERM_NUM 20000
#define BIG_NUM 10000

typedef struct {
	int c[MAX_DIM];
	double fit;
} Indiv;

typedef struct {
	int amount;
	Indiv v[4 * MAX_POP];
} Popu;

void randperm(int c[], int a, int b)
{
	int mark[MAX_DIM];
	int i, j, cnt, r;

	if (c == NULL || a > b) {
		printf("randperm error\n");
		exit(1);
	}

	for (i = a; i <= b; i++) {
		mark[i] = 0;
	}

	for (i = 0; i < b - a + 1; i++) {
		r = rand() % (b - a + 1 - i);
		cnt = 0;
		for (j = a; j < b + 1 && cnt < r - 1; j++) {
			if (mark[j] != 1) {
				cnt++;
			}
		}
		while (j < b + 1 && mark[j] == 1) {j++;}
		c[i] = j;
		mark[j] = 1;
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
				i++;
			}
			while (i < j && v[i].fit > tmp.fit) {
				i++;
			}
			if (i < j) {
				memcpy(&v[j], &v[i], sizeof(Indiv));
				j--;
			}
		}
		memcpy(&v[i], &tmp, sizeof(Indiv));
		fitSort(v, low, i - 1);
		fitSort(v, i + 1, high);
	}
}

double getFitness(const Mtx *m, int c[])
{
	int i, j;
	double fit;

	fit = 0;
	for (i = 0; i < m->dim - 1; i++) {
		fit += m->m[c[i]][c[i + 1]];
	}
	fit += m->m[c[m->dim - 1]][c[0]];

	return (1.0 / sqrt(fit));
}

void initPopu(const Mtx *m, Popu *p, int popSize)
{
	int i;

	for (i = 0; i < popSize; i++) {
		randperm(p->v[i].c, 0, m->dim - 1);
		p->v[i].fit = getFitness(m, p->v[i].c);
	}

	p->amount = popSize;
}

void speSelect(const Mtx *m, Popu *p, int maxPopSize)
{
	int i, j;
	int mark[MAX_POP * 4];
	double b;

	//calculate children's fitness
	for (i = maxPopSize; i < p->amount; i++) {
		p->v[i].fit = getFitness(m, p->v[i].c);
	}
	
	//sort fitness
	fitSort(p->v, 0, p->amount - 1);	

	//ingore repeat
	memset(mark, 0, sizeof(mark));
	for (i = 1; i < p->amount; i++) {
		if (p->v[i].fit == p->v[i - 1].fit) {
			mark[i] = 1;
		}
	}
	j = 1;
	for (i = 1; i < p->amount; i++) {
		if (mark[i] == 0 && j != i) {
			memcpy(&p->v[j], &p->v[i], sizeof(Indiv));
			j++;
		}
	}

	//select
	p->amount = maxPopSize;
}

void solveConflict(const Indiv *m, const Indiv *d, Indiv *child, int dim, int a, int b)
{
	int mark[MAX_DIM];
	int i, j, t;

	//memset(mark, 0, sizeof(mark));
	for (i = 0; i < dim; i++) {
		mark[i] = -1;
	}
	for (i = a; i<= b; i++) {
		mark[child->c[i]] = i;
	}
	for (i = 0; i < a; i++) {
		t = child->c[i];
		while (mark[t] != -1) {
			t = m->c[mark[t]];
		}
		child->c[i] = t;
	}
	for (i = b + 1; i < dim; i++) {
		t = child->c[i];
		while (mark[t] != -1) {
			t = m->c[mark[t]];
		}
		child->c[i] = t;
	}
}

void genChildren(const Indiv *m, const Indiv *d, Indiv *c1, Indiv *c2, int dim)
{
	int a, b, t;
	int i, j;

	a = rand() % dim;
	do {
		b = rand() % dim;
	} while (a == b);

	if (a > b) {
		t = a;
		a = b;
		b = t;
	}

	for (i = 0; i < a; i++) {
		c1->c[i] = m->c[i];
		c2->c[i] = d->c[i];
	}
	for (; i <= b; i++) {
		c1->c[i] = d->c[i];
		c2->c[i] = m->c[i];
	}
	for (; i < dim; i++) {
		c1->c[i] = m->c[i];
		c2->c[i] = d->c[i];
	}

	//solve conflict
	solveConflict(m, d, c1, dim, a, b);
	solveConflict(d, m, c2, dim, a, b);
}

void crossOver(const Mtx *m, Popu *p, double crossProb, double matingProb)
{
	int i, j, r;
	static int s[BIG_NUM * 2];
	double total;
	int a, b, childNum;
	int p1, p2;
	int e1, e2;

	total = 0;
	for (i = 0; i < p->amount; i++) {
		total += p->v[i].fit;
	}

	//calculate select prob
	a = 0;
	for (i = 0; i < p->amount; i++) {
		b = p->v[i].fit / total * (double)BIG_NUM + a;
		for (j = a; j < b; j++) {
			s[j] = i;
		}
		a = b;
	}

	//select mom and dad
	childNum = 0;
	for (i = 0; i < p->amount * crossProb; i++) {
		//select mom
		p1 = s[rand() % BIG_NUM];
		//select dad
		p2 = s[rand() % BIG_NUM];
		if (p1 == p2) {
			i--;
			continue;
		}
		//gen children
		if ((rand() % BIG_NUM / (double)BIG_NUM) < matingProb) {
			genChildren(&p->v[p1], &p->v[p2], &p->v[p->amount + childNum], \
							&p->v[p->amount + childNum + 1], m->dim);
			childNum += 2;
		}
	}
	
	p->amount += childNum;
}

void speMutat(const Mtx *m, Popu *p, int maxPopSize, double mutatProb)
{
	int i, j;
	int r1, r2, t;

	for (i = maxPopSize; i < p->amount; i++) {
		if ((rand() % BIG_NUM) / (double)BIG_NUM < mutatProb) {
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

	printf("Generation %d:\n", gen);
	printf("max fitness: %.10lf, min distance: %.6lf\n", p->v[0].fit, 1.0 / pow(p->v[0].fit, 2));
	printf("tour:\n");
	for (i = 0; i < m->dim; i++) {
		printf("%4d", p->v[0].c[i]);
		if ((i + 1) % 20 == 0) {
			printf("\n");
		}
	}
	printf("\n");
}

int main()
{ 
	Mtx mtx;
	Popu p;
	int i;

	allocMtx("att48.tsp", &mtx);
	if (mtx.m == NULL) {
		printf("alloc matrix error\n");
		return 1;
	}

	srand(time(NULL));

	initPopu(&mtx, &p, MAX_POP);
	for (i = 0; i < TERM_NUM; i++) {
		speSelect(&mtx, &p, MAX_POP);
		if ((i + 1) % 100 == 0) {
			showCurInfo(&mtx, &p, i);
		}
		crossOver(&mtx, &p, 0.9, 0.9);
		speMutat(&mtx, &p, MAX_POP, 0.2);
	}

	freeMtx(&mtx);
	return 0;
}
