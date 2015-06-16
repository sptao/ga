#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "tsp.h"
#include "parser.h"

#define MAX_DIM 300
#define POP_SIZE 256
#define TERM_NUM 1200000
#define BIG_NUM 10000
#define POP_NUM 4
#define MAX_PROC_NUM 256
#define MAX_FILE_LEN 4096

typedef struct {
	int c[MAX_DIM];
	double fit;
} Indiv;

typedef struct {
	int amount;
	Indiv v[4 * POP_SIZE];
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

int indivEqual(const Indiv *i1, const Indiv *i2, int dim)
{
	int i;
	int e;

	e = 1;
	for (i = 0; i < dim; i++) {
		if (i1->c[i] != i2->c[i]) {
			e = 0;
			break;
		}
	}

	return e;
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

	return (1.0 / fit);
}

void initPopu(const Mtx *m, Popu *p, int cnt)
{
	int i;

	for (i = 0; i < cnt; i++) {
		randperm(p->v[i].c, 0, m->dim - 1);
		p->v[i].fit = getFitness(m, p->v[i].c);
	}

	p->amount = cnt;
}

void speSelect(const Mtx *m, Popu *p, int maxPopSize)
{
	int i, j;
	int mark[POP_SIZE * 4];
	double b;

	//calculate children's fitness
//#pragma omp parallel for
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

void genChild(const Indiv *m, const Indiv *d, Indiv *c, int dim)
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
		c->c[i] = m->c[i];
	}
	for (; i <= b; i++) {
		c->c[i] = d->c[i];
	}
	for (; i < dim; i++) {
		c->c[i] = m->c[i];
	}

	//solve conflict
	solveConflict(m, d, c, dim, a, b);
}

void genChild1(const Indiv *m, const Indiv *d, Indiv *c, int dim)
{
	int a;
	int i, j;
	int mark[MAX_DIM];

	a = rand() % dim;
	memset(mark, 0, sizeof(mark));
	for (i = 0; i < a; i++) {
		c->c[i] = d->c[i]; 
		mark[d->c[i]] = 1;
	}
	for (j = 0; j < dim; j++) {
		if (mark[m->c[j]] != 1) {
			c->c[i++] = m->c[j];
		}
	}
}

void exterCrossOver(const Mtx *m, Popu *p1, Popu *p2, double crossProb, double matingProb)
{
	int i, j, r;
	static int s1[BIG_NUM * 2], s2[BIG_NUM * 2];
	double total1, total2;
	int a1, a2, b1, b2, childNum;
	int m1, m2;
	int e1, e2;

	total1 = 0;
	total2 = 0;
	for (i = 0; i < POP_SIZE; i++) {
		total1 += p1->v[i].fit;
		total2 += p2->v[i].fit;
	}

	//calculate select prob
	a1 = 0;
	a2 = 0;
	for (i = 0; i < POP_SIZE; i++) {
		b1 = p1->v[i].fit / total1 * (double)BIG_NUM + a1;
		for (j = a1; j < b1; j++) {
			s1[j] = i;
		}
		a1 = b1;
		b2 = p2->v[i].fit / total2 * (double)BIG_NUM + a2;
		for (j = a2; j < b2; j++) {
			s2[j] = i;
		}
		a2 = b2;
	}

	//select mom and dad
	j = p1->amount;
	for (i = 0; i < POP_SIZE * crossProb; i++) {
		//select mom
		m1 = s1[rand() % BIG_NUM];
		//select dad
		m2 = s2[rand() % BIG_NUM];
		//gen children
		if ((rand() % BIG_NUM / (double)BIG_NUM) < matingProb) {
			genChild1(&p1->v[m1], &p2->v[m2], &p1->v[j], m->dim);
			if (indivEqual(&p1->v[m1], &p1->v[j], m->dim) == 0) {
				j++;
			}
		}
	}
	
	p1->amount = j;
}

void interCrossOver(const Mtx *m, Popu *p, double crossProb, double matingProb)
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
	j = p->amount;
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
			genChild1(&p->v[p1], &p->v[p2], &p->v[j], m->dim);
			if (indivEqual(&p->v[p1], &p->v[j], m->dim) == 0 && \
							indivEqual(&p->v[p2], &p->v[j], m->dim) == 0) {
				j++;
			}
			genChild1(&p->v[p2], &p->v[p1], &p->v[j], m->dim);
			if (indivEqual(&p->v[p1], &p->v[j], m->dim) == 0 && \
							indivEqual(&p->v[p2], &p->v[j], m->dim) == 0) {
				j++;
			}
		}
	}
	
	p->amount = j;
}

void speMutat(const Mtx *m, Popu *p, int maxPopSize, double mutatProb)
{
	int i, j;
	int r1, r2, t;

	for (i = maxPopSize; i < p->amount; i++) {
		if ((rand() % BIG_NUM) / (double)BIG_NUM < mutatProb) {
			j = rand() % 4;
			while (j >= 0) {
				r1 = rand() % m->dim;
				r2 = rand() % m->dim;
				if (r1 != r2) {
					t = p->v[i].c[r1];
					p->v[i].c[r1] = p->v[i].c[r2];
					p->v[i].c[r2] = t;
				}

				j--;
			}
		}
	}
}

void speMutat1(const Mtx *m, Popu *p, int maxPopSize, double mutatProb)
{
	int i, j;
	int r1, r2, t;

	for (i = maxPopSize; i < p->amount; i++) {
		if ((rand() % BIG_NUM) / (double)BIG_NUM < mutatProb) {
			r1 = rand() % (m->dim - 1) ;
			t = p->v[i].c[r1];
			p->v[i].c[r1] = p->v[i].c[r1 + 1];
			p->v[i].c[r1 + 1] = t;
		}
	}
}

void showCurStatus(const Mtx *m, Popu *p, int gen, int pnum)
{
	int i;
	int cnt;

	printf("pnum: %d, Generation %d:\n", pnum, gen);
	printf("max fitness: %.10lf, min distance: %.6lf\n", p->v[0].fit, 1.0 / p->v[0].fit);

	cnt = 1;
	for (i = 1; i < p->amount; i++) {
		if (indivEqual(&p->v[i], &p->v[i - 1], m->dim) != 1) {
			cnt++;
		}
	}
	printf("different indiv count: %d\n", cnt);
	/*
	printf("tour:\n");
	for (i = 0; i < m->dim; i++) {
		printf("%4d", p->v[0].c[i]);
		if ((i + 1) % 20 == 0) {
			printf("\n");
		}
	}
	*/
	printf("\n");
}

int isMaster(int myid, int used_proc)
{
	if (myid % (used_proc / POP_NUM) == 0) {
		return 1;
	}

	return 0;
}

void genAlg(const Mtx *mtx, const MPI_Comm *cmm, int used_proc)
{
	Popu p;
	int ipp;
	int i, j, k;
	double h;

	ipp = POP_SIZE / (used_proc / POP_NUM);
	initPopu(mtx, &p, ipp);
	h = 0.000001;
	for (i = 0; i < TERM_NUM; i++) {
		speSelect(mtx, &p[j], POP_SIZE);
		intraCrossOver(mtx, &p[j], 0.9, 0.9);
		speMutat(mtx, &p[j], POP_SIZE, 0.3);

		if ((i + 1) % 2000 == 0) {
			if (isMaster()) {
				//get another pop's chromosome
			}
			//interCrossOver
			
			for (j = 0; j < POP_NUM; j++) {
				interCrossOver(mtx, &p[j], &p[(j + 1) % POP_NUM], 0.9, 0.9);
				showCurStatus(mtx, &p[j], i, j);
			}
			if (fabs(h - p[0].v[0].fit) / h < 0.01) {
				initPopu(mtx, &p[0], POP_SIZE);
			}
			else {
				h = p[0].v[0].fit;
			}
		}
	}
}

int main(int argc, char **argv)
{ 
	Mtx mtx;
	char buffer[MAX_FILE_LEN];
	int myid;
	int ipp, proc_num, used_proc, mutat_prob;
	MPI_Comm cmm;
	int i, color, key;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

	if (proc_num < POP_NUM) {
		printf("error: processor number should be no less than %d.", POP_NUM);
		exit(1);
	}
	used_proc = proc_num - proc_num % POP_NUM;

	//only process 0 read source tsp file
	if (myid == 0) {
		readFile("a280.tsp", buffer, MAX_FILE_LEN);
	}
	MPI_Bcast(buffer, MAX_FILE_LEN, MPI_CHAR, myid, MPI_COMM_WORLD);
	//each process generate the distance matrix
	allocMtx(buffer, &mtx);
	if (mtx.m == NULL) {
		printf("alloc matrix error\n");
		eixt(1);
	}

	srand(time(NULL));
	
	if (myid < used_proc) {
		color = i / (used_proc / POP_NUM);
		key = i;
	}
	else {
		color = MPI_UNDEFINED;
	}
	MPI_Comm_spilt(MPI_COMM_WORLD, color, key, &cmm);
	genAlg(&mtx, &cmm, used_proc);

	freeMtx(&mtx);
	MPI_Finilize();
	return 0;
}
