/*====================== display.c ========================
Contains functions for basic manipulation of a screen 
represented as a 2 dimensional array of colors.

A color is an ordered triple of ints, with each value standing
for red, green and blue respectively
==================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "gmath.h"
#include "ml6.h"
#include "display.h"
#include "matrix.h"

color change_color( int i ) {
  
  color c;
  i = i % 7;

  switch( i ) {
    
  case 0:
    c.red = 255;
    c.green = 255; 
    c.blue = 255;
    break;
  case 1:
    c.red = 255;
    c.green = 0; 
    c.blue = 0;
    break;
  case 2:
    c.red = 0;
    c.green = 255; 
    c.blue = 0;
    break;
  case 3:
    c.red = 0;
    c.green = 0; 
    c.blue = 255;
    break;
  case 4:
    c.red = 255;
    c.green = 255; 
    c.blue = 0;
    break;
  case 5:
    c.red = 255;
    c.green = 0; 
    c.blue = 255;
    break;
  case 6:
    c.red = 0;
    c.green = 255; 
    c.blue = 255;
    break;
  }
  return c;
}

/*======== void plot() ==========
Inputs:   screen s
         color c
         int x
         int y 
Returns: 
Sets the color at pixel x, y to the color represented by c
Note that s[0][0] will be the upper left hand corner 
of the screen. 
If you wish to change this behavior, you can change the indicies
of s that get set. For example, using s[x][YRES-1-y] will have
pixel 0, 0 located at the lower left corner of the screen

02/12/10 09:09:00
jdyrlandweaver
====================*/
void plot( screen s, color c, int x, int y, int z, struct matrix *zbuf) {
  int newy = YRES - 1 - y;
  if ( x >= 0 && x < XRES && newy >=0 && newy < YRES && z > (int)zbuf->m[x][newy]){
    s[x][newy] = c;
    zbuf->m[x][newy]=z;
  }
  //printf("plot?\n");
  // printf("plot color: %f, %f, %f\n", c.red, c.green, c.blue);
}
void plot1( screen s, int x, int y, int z, struct matrix *zbuf, double * n, struct constants *rcolor, color ambient, struct light ** point) {
  int newy = YRES - 1 - y;
  int i,j;
  color Ia,Id,Is;
  double *light_v = (double *)malloc(3*sizeof(double));
  double * view = (double *)malloc(3*sizeof(double));
  double * reflect = (double *)malloc(3*sizeof(double));
  double * normal = (double *)malloc(3*sizeof(double));
  normal[0]=n[0];
  normal[1]=n[1];
  normal[2]=n[2];
  view[0]=0;
  view[1]=0;
  view[2]=-1;
  color c;
  
  c.red = 0;
  c.green = 0;
  c.blue = 0;
  if ( x >= 0 && x < XRES && newy >=0 && newy < YRES && z > (int)zbuf->m[x][newy]){
    //printf("%f %f %f\n",n[0],n[1],n[2]);
    for(i = 0; point[i];i++){
      //ambient
      Ia.red = ambient.red * rcolor->r[0];
      Ia.green = ambient.green *rcolor->r[1];
      Ia.blue = ambient.blue *rcolor->r[2];
      //printf("%f %f %f\n",Ia.red,Ia.green,Ia.blue);
      //light vector
      light_v[0] = x - point[i]->l[0];
      light_v[1] = y - point[i]->l[1];
      light_v[2] = z - point[i]->l[2];

      //normalize
      normalize(light_v);
      normalize(normal);
      //printf("%f %f %f\n",n[0],n[1],n[2]);
      //diffuse
      Id.red = point[i]->c[0] * rcolor->g[0]*calculate_dot(light_v,normal)*-1;
      Id.green = point[i]->c[1] * rcolor->g[1]*calculate_dot(light_v,normal)*-1;
      Id.blue = point[i]->c[2] * rcolor->g[2]*calculate_dot(light_v,normal)*-1;
      //printf("light_v: %f %f %f\n",light_v[0],light_v[1],light_v[2]);
      //printf("normal: %f %f %f\n",n[0],n[1],n[2]);
      
      //specular
      reflect[0] = 2 * calculate_dot(light_v, normal)*normal[0] - light_v[0];
      reflect[1] = 2 * calculate_dot(light_v, normal)*normal[1] - light_v[1];
      reflect[2] = 2 * calculate_dot(light_v, normal)*normal[2] - light_v[2];
      //printf("%f\n",calculate_dot(light_v, n));
      //printf("%f %f %f\n",n[0],n[1],n[2]);
      Is.red = point[i]->c[0] * rcolor->b[0] * calculate_dot(reflect,view) * calculate_dot(reflect,view) * calculate_dot(reflect,view);
      Is.green = point[i]->c[1] * rcolor->b[1] * calculate_dot(reflect,view) * calculate_dot(reflect,view) * calculate_dot(reflect,view);
      Is.blue = point[i]->c[2] * rcolor->b[2] * calculate_dot(reflect,view) * calculate_dot(reflect,view) * calculate_dot(reflect,view);
      //printf("%f %f %f\n",reflect[0],reflect[1],reflect[2]);
      c.red += Ia.red + Id.red + Is.red;
      c.blue += Ia.blue + Id.blue + Is.blue;
      c.green += Ia.green + Id.green + Is.green;
      
      //printf("%f %f %f\n",c.red, c.green, c.blue);
      //printf("lighting\n");
      //printf("%d\n",i);
      /*c.red = c.red>255?255:c.red;
      c.red = c.red<0?0:c.red;
      c.green = c.green>255?255:c.green;
      c.green = c.green<0?0:c.green;
      c.blue = c.blue>255?255:c.blue;
      c.blue = c.blue<0?0:c.blue;*/
    }
    c.red = c.red>255?255:c.red;
      c.red = c.red<0?0:c.red;
      c.green = c.green>255?255:c.green;
      c.green = c.green<0?0:c.green;
      c.blue = c.blue>255?255:c.blue;
      c.blue = c.blue<0?0:c.blue;
      
    s[x][newy]=c;
    zbuf->m[x][newy]=z;
  }
  //printf("plot color: %f, %f, %f\n", c.red, c.green, c.blue);
}

/*======== void clear_screen() ==========
Inputs:   screen s  
Returns: 
Sets every color in screen s to black

02/12/10 09:13:40
jdyrlandweaver
====================*/
void clear_screen( screen s ) {

  int x, y;
  color c;

  c.red = 0;
  c.green = 0;
  c.blue = 0;

  for ( y=0; y < YRES; y++ )
    for ( x=0; x < XRES; x++)      
      s[x][y] = c;
}

/*======== void save_ppm() ==========
Inputs:   screen s
         char *file 
Returns: 
Saves screen s as a valid ppm file using the
settings in ml6.h

02/12/10 09:14:07
jdyrlandweaver
====================*/
void save_ppm( screen s, char *file) {

  int x, y;
  FILE *f;
  
  f = fopen(file, "w");
  fprintf(f, "P3\n%d %d\n%d\n", XRES, YRES, MAX_COLOR);
  for ( y=0; y < YRES; y++ ) {
    for ( x=0; x < XRES; x++) 
      
      fprintf(f, "%d %d %d ", (int) s[x][y].red, (int) s[x][y].green, (int) s[x][y].blue);
    fprintf(f, "\n");
  }
  fclose(f);
}
 
/*======== void save_extension() ==========
Inputs:   screen s
         char *file 
Returns: 
Saves the screen stored in s to the filename represented
by file. 
If the extension for file is an image format supported
by the "convert" command, the image will be saved in
that format.

02/12/10 09:14:46
jdyrlandweaver
====================*/
void save_extension( screen s, char *file) {
  
  int x, y;
  FILE *f;
  char line[256];

  sprintf(line, "convert - %s", file);

  f = popen(line, "w");
  fprintf(f, "P3\n%d %d\n%d\n", XRES, YRES, MAX_COLOR);
  for ( y=0; y < YRES; y++ ) {
    for ( x=0; x < XRES; x++)
      // if ((int) s[x][y].blue < 0 || (int) s[x][y].blue > 255)
      //   printf("color: %d, %d, %d\n",(int) s[x][y].red, (int) s[x][y].green, (int) s[x][y].blue);
      fprintf(f, "%d %d %d ", (int) s[x][y].red, (int) s[x][y].green, (int) s[x][y].blue);
    fprintf(f, "\n");
  }
  pclose(f);
}


/*======== void display() ==========
Inputs:   screen s 
Returns: 
Will display the screen s on your monitor

02/12/10 09:16:30
jdyrlandweaver
====================*/
void display( screen s) {
   int x, i;
  char *fname = ".tmp.png";
  save_extension(s, fname);
  i = fork();
  if (i == 0) {
    execlp("display", "display", fname, NULL);
  }
  else {
    wait(&x);
    remove( fname );
  }
  /* For some reason, this refuses to run correctly
     on some systems. Most likely a strange imagemagick
     install issue. 
     Above is a workaroudn for now.
  int x, y;
  FILE *f;

  f = popen("display", "w");

  fprintf(f, "P3\n%d %d\n%d\n", XRES, YRES, MAX_COLOR);
  for ( y=0; y < YRES; y++ ) {
    for ( x=0; x < XRES; x++) 
      
      fprintf(f, "%d %d %d ", s[x][y].red, s[x][y].green, s[x][y].blue);
    fprintf(f, "\n");
  }
  pclose(f);
  */
}

