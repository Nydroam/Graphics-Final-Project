#ifndef GMATH_H
#define GMATH_H

#include "matrix.h"

double * calculate_normal( double a1, double a2, double a3,
			   double b1, double b2, double b3 );
//double calculate_dot( struct matrix *points, int i );
double calculate_dot(double * v1, double * v2);

void * normalize(double * v);

int nearly_equal(double n1, double n2);
#endif
