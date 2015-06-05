#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tsp.h"

const char *strType = "TYPE";
const char *strDim = "DIMENSION";
const char *strE2d = "EUC_2D";
const char *strSec = "NODE_COORD_SECTION";

void allocMtx(const char *filename, Mtx *mtx)
{
	if (mtx == NULL) {
		return;
	}

	mtx->m = NULL;
	FILE *fp = fopen(filename, "rb");
	if (fp == NULL) {
		printf("open file error\n");
		return;
	}
	
	int len;
	char *buf, *ptr;
	double *x, *y;
	int i, j;

	fseek(fp, 0, SEEK_END);
	len = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	buf = (char *)malloc(sizeof(char) * len);
	if (buf == NULL) {
		fclose(fp);
		return;
	}
	len = fread(buf, 1, len, fp);
	fclose(fp);

	//need to do file format check
	
	ptr = strstr(buf, strDim);
	if (ptr == NULL) {
		free(buf);
		printf("incorrect format\n");
		return;
	}
	while (!(*ptr >= '0' && *ptr <= '9')) {ptr++;}	//find the dimension
	sscanf(ptr, "%d", &mtx->dim);
	ptr = strstr(ptr, strSec);						//go to line: NODE_COORD_SECTION
	if (ptr == NULL) {
		free(buf);
		printf("incorrect format\n");
		return;
	}
	while (*ptr != '\n') {ptr++;}	
	ptr++;			//go to the data line
	x = (double *)malloc(sizeof(double) * mtx->dim);
	y = (double *)malloc(sizeof(double) * mtx->dim);
	if (x == NULL || y == NULL) {
		printf("alloc memory error\n");
		if (!x) {
			free(x);
		}
		if (!y) {
			free(x);
		}
		free(buf);
		return;
	}
	for (i = 0; i < mtx->dim; i++) {
		sscanf(ptr, "%d%lf%lf", &j, &x[i], &y[i]);
		while (*ptr != '\n') {ptr++;}
		ptr++;		//go to the next line
	}

	mtx->m = (double **)malloc(sizeof(double *) * mtx->dim);
	if (mtx->m == NULL) {
		printf("alloc memory error\n");
		free(buf);
		return;
	}
	mtx->m[0] = (double *)malloc(sizeof(double) * mtx->dim * mtx->dim);
	if (mtx->m[0] == NULL) {
		printf("alloc memory error\n");
		free(mtx->m);
		free(buf);
		mtx->m = NULL;
		return;
	}

	for (i = 0; i < mtx->dim; i++) {
		mtx->m[i] = mtx->m[0] + i * mtx->dim;
		for (j = 0; j < i; j++) {
			mtx->m[i][j] = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));
		}
		mtx->m[i][i] = 0;
	}
	for (i = 0; i < mtx->dim; i++) {
		for (j = i + 1; j < mtx->dim; j++) {
			mtx->m[i][j] = mtx->m[j][i];
		}
	}

	free(x);
	free(y);
	free(buf);
	return;
}

void freeMtx(Mtx *mtx)
{
	if (mtx != NULL && mtx->m != NULL) {
		free(mtx->m[0]);
		free(mtx->m);
	}
}
