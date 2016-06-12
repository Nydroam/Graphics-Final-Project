#ifndef DISPLAY_H
#define DISPLAY_H

#include "symtab.h"
#include "matrix.h"
#include "gmath.h"
void plot( screen s, color c, int x, int y, int z, struct matrix* zbuf);
void plot1( screen s, int x, int y, int z, struct matrix* zbuf, double *n, struct constants * rcolor, color ambient, struct light ** point);
void clear_screen( screen s);
void save_ppm( screen s, char *file);
void save_extension( screen s, char *file);
void display( screen s);

color change_color( int i );
#endif
