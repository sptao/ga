#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "tsp.h"
#include "parser.h"

typedef char chromoType;

#define MAX_DIM 256
#define POP_SIZE 256
#define TERM_NUM 1200000
#define BIG_NUM 10000
#define POP_NUM 4
#define MAX_PROC_NUM 256
#define MAX_FILE_LEN 4096

typedef struct {
	chromoType c[POP_SIZE][MAX_DIM];
	double fit[POP_SIZE];
	int amount;
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

int indivEqual(const int c1[], const int c2[], int dim)
{
	int i;
	int e;

	e = 1;
	for (i = 0; i < dim; i++) {
		if (c1[i] != c2[i]) {
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

void initPopu(const Mtx *m, Popu *p, int ipp)
{
	int i;

	for (i = 0; i < ipp; i++) {
		randperm(p->c[i], 0, m->dim - 1);
		p->fit[i] = getFitness(m, p->c[i]);
	}

	p->amount = ipp;
}

void speSelect(const Mtx *m, Popu *p, int used_proc, MPI_Comm cmm)
{
	Popu ap;
	int i, j, npp, ipp, rank;
	int cnt[POP_SIZE];
	int disp[POP_SIZE];
	double fit[POP_SIZE], mid;
	double b;

	//number of processes per population
	npp = used_proc / POP_NUM;
	//individual number per process
	ipp = POP_SIZE / npp;
	//get process's rank in each population
	MPI_Comm_rank(cmm, &rank);

	//calculate children's fitness
	for (i = ipp; i < p->amount; i++) {
		p->fit[i] = getFitness(m, p->c[i]);
	}
	
	//send fit to master
	MPI_Gather(p->fit, p->amount, MPI_DOUBLE, fit, \
					POP_SIZE, MPI_DOUBLE, 0, cmm);
	//master process find the mid value and broadcast to all processes
	if (rank == 0) {
		mid = findMidFit(fit, );
	}
	MPI_Bcast(&mid, 1, MPI_DOUBLE, 0, cmm);

	ap.amount = 0;
	for (i = 0; i < p->amount; i++) {
		if (p->fit[i] >= mid) {
			memcpy(ap.c[ap.amount], p->c[i], sizeof(chromoType) * MAX_DIM);
			ap.fit[ap.amount] = p->fit[i];
			ap.amount++;
		}
	}

	//gather alive count 
	MPI_Allgather(&ap.amount, 1, MPI_INT, cnt, npp, MPI_INT, cmm);
	//gather alive fit 
	disp[0] = 0;
	for (i = 1; i < npp; i++) {
		disp[i] = disp[i - 1] + cnt[i - 1];
	}
	MPI_Allgatherv(ap.fit, ap.amount, MPI_DOUBLE, p->fit, cnt, disp, MPI_DOUBLE, cmm);
	//gather alive chromosome
	disp[0] = 0;
	for (i = 1; i < npp; i++) {
		disp[i] = disp[i - 1] + cnt[i - 1] * sizeof(chromoType) * MAX_DIM;
	}
	MPI_Allgatherv(ap.c[0], ap.amount * sizeof(chromoType) * MAX_DIM, \
					MPI_CHAR, p->c[0], cnt * sizeof(chromoType) * MAX_DIM, disp, MPI_CHAR, cmm);

	//set popu amount
	p->amount = POP_SIZE;
}

void solveConflict(const int m[], const int d[], int child[], int dim, int a, int b)
{
	int mark[MAX_DIM];
	int i, j, t;

	//memset(mark, 0, sizeof(mark));
	for (i = 0; i < dim; i++) {
		mark[i] = -1;
	}
	for (i = a; i<= b; i++) {
		mark[child[i]] = i;
	}
	for (i = 0; i < a; i++) {
		t = child[i];
		while (mark[t] != -1) {
			t = m[mark[t]];
		}
		child[i] = t;
	}
	for (i = b + 1; i < dim; i++) {
		t = child[i];
		while (mark[t] != -1) {
			t = m[mark[t]];
		}
		child[i] = t;
	}
}

void genChild1(const int m[], const int d[], int c[], int dim)
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
		c[i] = m[i];
	}
	for (; i <= b; i++) {
		c[i] = d[i];
	}
	for (; i < dim; i++) {
		c[i] = m[i];
	}

	//solve conflict
	solveConflict(m, d, c, dim, a, b);
}

void genChild1(const int m[], const int d[], int c[], int dim)
{
	int a;
	int i, j;
	int mark[MAX_DIM];

	a = rand() % dim;
	memset(mark, 0, sizeof(mark));
	for (i = 0; i < a; i++) {
		c[i] = d[i]; 
		mark[d[i]] = 1;
	}
	for (j = 0; j < dim; j++) {
		if (mark[m[j]] != 1) {
			c[i++] = m[j];
		}
	}
}

void interCrossOver(const Mtx *m, Popu *p1, Popu *p2, int ipp, MPI_Comm cmm)
{
	int i, j, r;
	static int s1[BIG_NUM * 2], s2[BIG_NUM * 2];
	Popu ap;
	int rank;
	double total1, total2;
	int a1, a2, b1, b2, childNum;
	int m1, m2;

	total1 = 0;
	total2 = 0;
	for (i = 0; i < MAX_POP; i++) {
		total1 += p1->v[i].fit;
		total2 += p2->v[i].fit;
	}

	//calculate select prob
	a1 = 0;
	a2 = 0;
	for (i = 0; i < MAX_POP; i++) {
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

	//each processor produce same amount of children
	ap.amount = 0;
	while (ap.amount < ipp) {
		m1 = s1[rand() % BIG_NUM];
		m2 = s2[rand() % BIG_NUM];
		if (m1 == m2) {
			continue;
		}
		genChild1(p1->c[m1], p2->c[m2], ap.c[ap.amount], m->dim);
		if (indivEqual(p1->c[m1], ap.c[ap.amount], m->dim) == 0 && \
						indivEqual(p2->c[m2], ap.c[ap.amount], m->dim) == 0) {
			ap.amount++;
		}
	}
	MPI_Comm_rank(cmm, &rank);
	memcpy(p->c[0], p->c[rank * ipp], sizeof(chromoType) * MAX_DIM * ipp);
	memcpy(p->fit, &p->fit[rank * ipp], sizeof(double) * ipp);
	memcpy(p->c[ipp], ap.c[0], sizeof(chromoType) * MAX_DIM * ipp);
	memcpy(&p->fit[ipp], ap.fit, sizeof(double) * ipp);
}

void intraCrossOver(const Mtx *m, Popu *p, int ipp, MPI_Comm cmm)
{
	int i, j;
	static int s[BIG_NUM * 2];
	Popu ap;
	int rank;
	double total;
	int a, b;
	int m1, m2;

	total = 0;
	for (i = 0; i < p->amount; i++) {
		total += p->fit[i];
	}

	//calculate select prob
	a = 0;
	for (i = 0; i < p->amount; i++) {
		b = p->fit[i] / total * (double)BIG_NUM + a;
		for (j = a; j < b; j++) {
			s[j] = i;
		}
		a = b;
	}

	//each processor produce same amount of children
	for (ap.amount = 0; ap.amount < ipp;) {
		m1 = s[rand() % BIG_NUM];
		m2 = s[rand() % BIG_NUM];
		if (m1 == m2) {
			continue;
		}
		genChild1(p->c[m1], p->c[m2], ap.c[ap.amount], m->dim);
		if (indivEqual(p->c[m1], ap.c[ap.amount], m->dim) == 0 && \
						indivEqual(p->c[m2], ap.c[ap.amount], m->dim) == 0) {
			ap.amount++;
		}
	}
	MPI_Comm_rank(cmm, &rank);
	memcpy(p->c[0], p->c[rank * ipp], sizeof(chromoType) * MAX_DIM * ipp);
	memcpy(p->fit, &p->fit[rank * ipp], sizeof(double) * ipp);
	memcpy(p->c[ipp], ap.c[0], sizeof(chromoType) * MAX_DIM * ipp);
	memcpy(&p->fit[ipp], ap.fit, sizeof(double) * ipp);
}

void speMutat(const Mtx *m, Popu *p, int ipp, double mutatProb)
{
	int i, j;
	int r1, r2, t;

	for (i = p->amount - ipp; i < p->amount; i++) {
		if ((rand() % BIG_NUM) / (double)BIG_NUM < mutatProb) {
			j = rand() % 4;
			while (j >= 0) {
				r1 = rand() % m->dim;
				r2 = rand() % m->dim;
				if (r1 != r2) {
					t = p->c[r1][i];
					p->c[r1][i] = p->c[r2][i];
					p->c[r2][i] = t;
				}

				j--;
			}
		}
	}
}

void showCurStatus(const Mtx *m, Popu *p, int gen, int pnum)
{
	int i;
	int cnt;

	printf("pop id: %d, Generation %d:\n", pnum, gen);
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

void genAlg(const Mtx *mtx, const MPI_Comm *cmm, int myid, int used_proc)
{
	Popu p, ap;
	int ipp;
	int i, j, k, status;
	double h;

	ipp = POP_SIZE / (used_proc / POP_NUM);
	initPopu(mtx, &p, ipp);
	h = 0.000001;
	for (i = 0; i < TERM_NUM; i++) {
		speSelect(mtx, &p, used_proc, cmm)
		if ((i + 1) % 2000 != 0) {
			intraCrossOver(mtx, &p, ipp, cmm)
			speMutat(mtx, &p, ipp, 0.3)
		}
		else {
			if (myid / (used_proc / POP_NUM) != POP_NUM - 1) {
				MPI_Send(&p, sizeof(Popu), MPI_CHAR, myid + used_proc / POP_NUM, 0, MPI_COMM_WORLD);
				MPI_Recv(&ap, sizeof(Popu), MPI_CHAR, \
							(myid + used_proc - used_proc / POP_NUM) % used_proc, \
							0, MPI_COMM_WORLD, &status);
			}
			else {
				MPI_Recv(&ap, sizeof(Popu), MPI_CHAR, \
							(myid + used_proc - used_proc / POP_NUM) % used_proc, \
							0, MPI_COMM_WORLD, &status);
				MPI_Send(&p, sizeof(Popu), MPI_CHAR, myid + used_proc / POP_NUM, 0, MPI_COMM_WORLD);
			}
			interCrossOver(mtx, &p, &ap, ipp, cmm)
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
		color = myid / (used_proc / POP_NUM);
		key = myid;
	}
	else {
		color = MPI_UNDEFINED;
	}
	MPI_Comm_spilt(MPI_COMM_WORLD, color, key, &cmm);
	genAlg(mtx, cmm, myid, used_proc)

	freeMtx(&mtx);
	MPI_Finilize();
	return 0;
}
