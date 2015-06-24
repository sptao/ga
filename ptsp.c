#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "tsp.h"
#include "parser.h"

#define MAX_DIM 256
#define POP_SIZE 256
#define TERM_NUM 1200000
#define BIG_NUM 10000
#define POP_NUM 2
#define MAX_PROC_NUM 256
#define MAX_FILE_LEN 8192
#define INTER_CROSSOVER_NUM 2000

typedef int chromoType;
typedef struct {
	chromoType c[POP_SIZE * 2][MAX_DIM];
	double fit[POP_SIZE * 2];
	int amount;
} Popu;

void randperm(chromoType c[], int a, int b)
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

int indivEqual(const chromoType c1[], const chromoType c2[], int dim)
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

double getFitness(const Mtx *m, chromoType c[])
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

void sortFit(double *f, int low, int high)
{
	int i, j;
	double tmp;

	if (low < high) {
		i = low;
		j = high;
		tmp = f[i];
		while (i < j) {
			while (i < j && f[j] < tmp) {
				j--;
			}
			if (i < j) {
				f[i] = f[j];
				i++;
			} 
			while (i < j && f[i] > tmp) {
				i++;
			}
			if (i < j) {
				f[j] = f[i];
				j--;
			}
		}
		f[i] = tmp;
		sortFit(f, low, i - 1);
		sortFit(f, i + 1, high);
	}
}

double findKthFit(double *fit, int fitCnt, int k)
{
	sortFit(fit, 0, fitCnt - 1);
	
	return fit[k - 1];
}

void speSelect(const Mtx *m, Popu *p, int used_proc, MPI_Comm cmm)
{
	Popu ap;
	int i, j, npp, ipp, rank, fitCnt, eq;
	int cnt[POP_SIZE];
	int disp[POP_SIZE];
	double fit[POP_SIZE * 2], kth;
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
	//printf("start select\n");
	MPI_Reduce(&p->amount, &fitCnt, 1, MPI_INT, MPI_SUM, 0, cmm);
	//printf("rank = %d, p->amount = %d\n", rank, p->amount);
	MPI_Gather(p->fit, p->amount, MPI_DOUBLE, fit, \
					p->amount, MPI_DOUBLE, 0, cmm);
	//printf("finish gather fit\n");
	//master process find the kth fit value and broadcast to all processes
	if (rank == 0) {
		//printf("start find, fitCnt = %d\n", fitCnt);
		kth = findKthFit(fit, fitCnt, POP_SIZE);
		//printf("finish find\n");
	}
	//printf("start bcast: %d\n", rank);
	MPI_Bcast(&kth, 1, MPI_DOUBLE, 0, cmm);
	//printf("finish bcast\n");

	//store selected indiv in another Popu
	ap.amount = 0;
	eq = 0;
	for (i = 0; i < p->amount; i++) {
		if (p->fit[i] > kth) {
			memcpy(ap.c[ap.amount], p->c[i], sizeof(chromoType) * MAX_DIM);
			ap.fit[ap.amount] = p->fit[i];
			ap.amount++;
		}
		else if (p->fit[i] == kth && eq == 0) {
			memcpy(ap.c[ap.amount], p->c[i], sizeof(chromoType) * MAX_DIM);
			ap.fit[ap.amount] = p->fit[i];
			ap.amount++;
			eq = 1;
		}
	}

	//gather remainder count 
	//printf("start allgather amount\n");
	MPI_Allgather(&ap.amount, 1, MPI_INT, cnt, 1, MPI_INT, cmm);
	//set popu amount
	p->amount = 0; 
	for (i = 0; i < npp; i++) {
		p->amount += cnt[i];
	}
	//gather remainder fit 
	disp[0] = 0;
	for (i = 1; i < npp; i++) {
		disp[i] = disp[i - 1] + cnt[i - 1];
	}
	//printf("start gatherv fit\n");
	MPI_Allgatherv(ap.fit, ap.amount, MPI_DOUBLE, p->fit, cnt, disp, MPI_DOUBLE, cmm);
	//printf("gather fit over\n");
	//gather remainder chromosome
	disp[0] = 0;
	cnt[0] = cnt[0] * sizeof(chromoType) * MAX_DIM;
	for (i = 1; i < npp; i++) {
		cnt[i] = cnt[i] * sizeof(chromoType) * MAX_DIM;
		disp[i] = disp[i - 1] + cnt[i - 1];
	}
	//printf("start gatherv c\n");
	MPI_Allgatherv(&ap.c[0][0], ap.amount * sizeof(chromoType) * MAX_DIM, \
					MPI_CHAR, &p->c[0][0], cnt, disp, MPI_CHAR, cmm);		
	//printf("gatherv over\n");
}

void solveConflict(const chromoType m[], const chromoType d[], chromoType child[], int dim, int a, int b)
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

void genChild(const chromoType m[], const chromoType d[], chromoType c[], int dim)
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

void genChild1(const chromoType m[], const chromoType d[], chromoType c[], int dim)
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
	for (i = 0; i < POP_SIZE; i++) {
		total1 += p1->fit[i];
		total2 += p2->fit[i];
	}

	//calculate select prob
	a1 = 0;
	a2 = 0;
	for (i = 0; i < POP_SIZE; i++) {
		b1 = p1->fit[i] / total1 * (double)BIG_NUM + a1;
		for (j = a1; j < b1; j++) {
			s1[j] = i;
		}
		a1 = b1;
		b2 = p2->fit[i] / total2 * (double)BIG_NUM + a2;
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
	memcpy(p1->c[0], p1->c[rank * ipp], sizeof(chromoType) * MAX_DIM * ipp);
	memcpy(p1->fit, &p1->fit[rank * ipp], sizeof(double) * ipp);
	memcpy(p1->c[ipp], ap.c[0], sizeof(chromoType) * MAX_DIM * ipp);
	memcpy(&p1->fit[ipp], ap.fit, sizeof(double) * ipp);
}

void intraCrossOver(const Mtx *m, Popu *p, int ipp, MPI_Comm cmm)
{
	int i, j;
	static int s[BIG_NUM * 2];
	Popu ap;
	int rank;
	double total;
	int a;
	int m1, m2;

	total = 0;
	for (i = 0; i < p->amount; i++) {
		total += p->fit[i];
	}

	//calculate select prob
	a = 0;
	j = 0;
	for (i = 0; i < p->amount; i++) {
		a += p->fit[i] / total * (double)BIG_NUM;
		while (j < a) {
			s[j] = i;
			j++;
		}
	}

	//each processor produce same amount of children
	int cnt = 0;
	for (ap.amount = 0; ap.amount < ipp;) {
		m1 = s[rand() % a];
		m2 = s[rand() % a];
		if (m1 == m2) {
			continue;
		}
		genChild(p->c[m1], p->c[m2], ap.c[ap.amount], m->dim);
		if (indivEqual(p->c[m1], ap.c[ap.amount], m->dim) == 0 && \
						indivEqual(p->c[m2], ap.c[ap.amount], m->dim) == 0) {
			ap.amount++;
		}
		cnt++;
		if (cnt > 1000) {
			break;
		}
	}
	while (ap.amount < ipp) {
		randperm(ap.c[ap.amount], 0, m->dim - 1);
		ap.amount++;
	}

	MPI_Comm_rank(cmm, &rank);
	if (rank != 0) {
		memcpy(p->c[0], p->c[rank * ipp], sizeof(chromoType) * MAX_DIM * ipp);
		memcpy(p->fit, &p->fit[rank * ipp], sizeof(double) * ipp);
	}
	memcpy(p->c[ipp], ap.c[0], sizeof(chromoType) * MAX_DIM * ipp);
	memcpy(&p->fit[ipp], ap.fit, sizeof(double) * ipp);
	p->amount = 2 * ipp;
}

void speMutat(const Mtx *m, Popu *p, int ipp, double mutatProb)
{
	int i, j;
	int r1, r2, t;

	for (i = ipp; i < p->amount; i++) {
		if ((rand() % BIG_NUM) / (double)BIG_NUM < mutatProb) {
			j = rand() % 2 + 1;
			while (j >= 0) {
				r1 = rand() % m->dim;
				r2 = rand() % m->dim;
				if (r1 != r2) {
					t = p->c[i][r1];
					p->c[i][r1] = p->c[i][r2];
					p->c[i][r2] = t;
				}

				j--;
			}
		}
	}
}

void showCurStatus(const Mtx *m, Popu *p, int gen, int npp, int myid)
{
	int i, n;
	int cnt;

	if (myid % npp == 0) {
		printf("pop id: %d, Generation %d:\n", myid / npp, gen);
		n = 0;
		if (p->amount != 256) {
			printf("amount error: %d\n", p->amount);
		}
		for (i = 1; i < p->amount; i++) {
			if (p->fit[i] > p->fit[n]) {
				n = i;
			}
		}
		printf("max fitness: %.10lf, min distance: %.6lf\n", p->fit[n], 1.0 / p->fit[n]);
	}

	if (myid == 0) {
		printf("\n");
	}
}

void genAlg(const Mtx *mtx, MPI_Comm cmm, int used_proc)
{
	Popu p, ap;
	int myid, ipp, npp;
	int i, j, k;
	MPI_Status status;
	double h;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	npp = used_proc / POP_NUM;
	ipp = POP_SIZE / npp;
	printf("start init. ipp = %d\n", ipp);
	initPopu(mtx, &p, ipp);
	printf("init over\n");
	h = 0.000001;
	for (i = 0; i < TERM_NUM; i++) {
		//species select according to fitness
		speSelect(mtx, &p, used_proc, cmm);
		if ((i + 1) % INTER_CROSSOVER_NUM != 0) {
			//crossover in one population
			intraCrossOver(mtx, &p, ipp, cmm);
			//printf("crossover: %d\n", i);
		}
		else {
			//show current status
			showCurStatus(mtx, &p, i, npp, myid);
			//exchange info between Popus
			if (POP_NUM > 1) {
				if (myid / npp != POP_NUM - 1) {
					MPI_Send(&p, sizeof(Popu), MPI_CHAR, (myid + npp) % used_proc, 0, MPI_COMM_WORLD);
					MPI_Recv(&ap, sizeof(Popu), MPI_CHAR, \
								(myid + used_proc - npp) % used_proc, \
								0, MPI_COMM_WORLD, &status);
				}
				else {
					MPI_Recv(&ap, sizeof(Popu), MPI_CHAR, \
								(myid + used_proc - npp) % used_proc, \
								0, MPI_COMM_WORLD, &status);
					MPI_Send(&p, sizeof(Popu), MPI_CHAR, (myid + npp) % used_proc, 0, MPI_COMM_WORLD);
				}
				//crossover in two popultaion
				interCrossOver(mtx, &p, &ap, ipp, cmm);
			}
		}

		//mutation
		speMutat(mtx, &p, ipp, 0.7);
	}
}

int main(int argc, char **argv)
{ 
	Mtx mtx;
	char buffer[MAX_FILE_LEN];
	int myid;
	int npp, proc_num, used_proc;
	MPI_Comm cmm;
	int color, key;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

	//process numbers must be less than POP_NUM
	if (proc_num < POP_NUM) {
		if (myid == 0) {
			printf("error: processor number should be no less than %d.\n", POP_NUM);
		}
		MPI_Finalize();
		return 0;
	}
	//some processes do not work
	used_proc = proc_num - proc_num % POP_NUM;

	//only process 0 read source tsp file
	if (myid == 0) {
		readFile("eil51.tsp", buffer, MAX_FILE_LEN);
	}
	MPI_Bcast(buffer, MAX_FILE_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);
	//each process generate the same distance matrix
	allocMtx(buffer, &mtx);
	if (mtx.m == NULL) {
		printf("alloc matrix error\n");
		exit(1);
	}

	srand(time(NULL) + myid);
	
	//create new communicator
	npp = used_proc / POP_NUM;		//number of processes per population
	if (myid < used_proc) {
		color = myid / npp;
		key = myid % npp;
	}
	else {
		color = MPI_UNDEFINED;
	}
	MPI_Comm_split(MPI_COMM_WORLD, color, key, &cmm);

	//start algorithm
	printf("start gen alg\nmtx.dim = %d\n", mtx.dim);
	genAlg(&mtx, cmm, used_proc);

	freeMtx(&mtx);
	MPI_Finalize();
	return 0;
}
