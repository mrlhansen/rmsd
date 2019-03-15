#ifndef UTILS_H
#define UTILS_H

#define error(msg) \
	fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, msg); \
	exit(-1)

void find_dbl(int, char**, char*, double*);
void find_int(int, char**, char*, int*);
void find_str(int, char**, char*, char*);
void find_opt(int, char**, char*, int*);

void randu(double*, int);
void randi(int*, int, int, int);
void randg(double*, int);

#endif
