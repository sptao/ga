#ifndef __PARSER_H__
#define __PARSER_H__

#include "tsp.h"

void readFile(const char *filename, char *buf, int bufLen);
void allocMtx(const char *buf, Mtx *mtx);
void freeMtx(Mtx *mtx);

#endif
