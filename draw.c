#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"
#include "gmath.h"
#include "symtab.h"

/*======== void add_polygon() ==========
Inputs:   struct matrix *surfaces
				 double x0
				 double y0
				 double z0
				 double x1
				 double y1
				 double z1
				 double x2
				 double y2
				 double z2  
Returns: 
Adds the vertices (x0, y0, z0), (x1, y1, z1)
and (x2, y2, z2) to the polygon matrix. They
define a single triangle surface.

04/16/13 13:05:59
jdyrlandweaver
====================*/
void add_polygon( struct matrix *polygons, 
			double x0, double y0, double z0, 
			double x1, double y1, double z1, 
			double x2, double y2, double z2 ) {
	add_point(polygons, x0, y0, z0);
	add_point(polygons, x1, y1, z1);
	add_point(polygons, x2, y2, z2);
}
	

	/*
			struct constants rcolor;
			rcolor.r[LAMBIENT] = 1;
			rcolor.g[LAMBIENT] = 1;
			rcolor.b[LAMBIENT] = 1;
			
			rcolor.r[LDIFFUSE] = 0;
			rcolor.g[LDIFFUSE] = 0;
			rcolor.b[LDIFFUSE] = 0;

			rcolor.r[LSPECULAR] = 0;
			rcolor.g[LSPECULAR] = 0;
			rcolor.b[LSPECULAR] = 0;
						color ambient;
			ambient.red = 255;
			ambient.green = 255;
			ambient.blue = 255;
			*/

/*======== void draw_polygons() ==========
Inputs:   struct matrix *polygons
					screen s
					color c  
Returns: 
Goes through polygons 3 points at a time, drawing 
lines connecting each points to create bounding
triangles

04/16/13 13:13:27
jdyrlandweaver
====================*/
void draw_polygons( struct matrix * polygons, screen s, color c, struct matrix* zbuf, struct constants *rcolor, color ambient, struct light *point) {
	
	int i,x,y;  
	double xB, yB, xM, yM, xT, yT, d0, d1, d2;
	double ax, ay, az, bx, by, bz;
	double * normal;
	double * view;
	double * light_v;
	double * reflect;
	double xL,xR;
	double zB, zM, zT;
	double zL,zR;
	color Ia, Id, Is;
	double dot;
	view = (double *)malloc(3 * sizeof(double));
	light_v = (double *)malloc(3 * sizeof(double));
	reflect = (double *)malloc(3 * sizeof(double));
	view[0] = 0;
	view[1] = 0;
	view[2] = -1;
	for( i=0; i < polygons->lastcol-2; i+=3 ) {
		//get the surface normal
		ax = polygons->m[0][i+1] - polygons->m[0][i];
		ay = polygons->m[1][i+1] - polygons->m[1][i];
		az = polygons->m[2][i+1] - polygons->m[2][i];
		bx = polygons->m[0][i+2] - polygons->m[0][i];
		by = polygons->m[1][i+2] - polygons->m[1][i];
		bz = polygons->m[2][i+2] - polygons->m[2][i];
		normal = calculate_normal( ax, ay, az, bx, by, bz );

		if ( calculate_dot( normal, view ) < 0 ) {
			//determines top mid bot
			if (polygons->m[1][i] >= polygons->m[1][i+1] && polygons->m[1][i] >= polygons->m[1][i+2]) {
				xT = polygons->m[0][i];
				yT = polygons->m[1][i];
				zT = polygons->m[2][i];
				if (polygons->m[1][i+1] >= polygons->m[1][i+2]) {
					xM = polygons->m[0][i+1];
					yM = polygons->m[1][i+1];
					zM = polygons->m[2][i+1];
					xB = polygons->m[0][i+2];
					yB = polygons->m[1][i+2];
					zB = polygons->m[2][i+2];
				}
				else {
					xB = polygons->m[0][i+1];
					yB = polygons->m[1][i+1];
					zB = polygons->m[2][i+1];
					xM = polygons->m[0][i+2];
					yM = polygons->m[1][i+2];
					zM = polygons->m[2][i+2];
				}
			}
			else if (polygons->m[1][i+1] >= polygons->m[1][i] && polygons->m[1][i+1] >= polygons->m[1][i+2]) {
				xT = polygons->m[0][i+1];
				yT = polygons->m[1][i+1];
				zT = polygons->m[2][i+1];
				if (polygons->m[1][i] >= polygons->m[1][i+2]) {
					xM = polygons->m[0][i];
					yM = polygons->m[1][i];
					zM = polygons->m[2][i];
					xB = polygons->m[0][i+2];
					yB = polygons->m[1][i+2];
					zB = polygons->m[2][i+2];
				}
				else {
					xB = polygons->m[0][i];
					yB = polygons->m[1][i];
					zB = polygons->m[2][i];
					xM = polygons->m[0][i+2];
					yM = polygons->m[1][i+2];
					zM = polygons->m[2][i+2];
				}
			}
			else {
				xT = polygons->m[0][i+2];
				yT = polygons->m[1][i+2];
				zT = polygons->m[2][i+2];
				if (polygons->m[1][i] >= polygons->m[1][i+1]) {
					xM = polygons->m[0][i];
					yM = polygons->m[1][i];
					zM = polygons->m[2][i];
					xB = polygons->m[0][i+1];
					yB = polygons->m[1][i+1];
					zB = polygons->m[2][i+1];
				}
				else {
					xB = polygons->m[0][i];
					yB = polygons->m[1][i];
					zB = polygons->m[2][i];
					xM = polygons->m[0][i+1];
					yM = polygons->m[1][i+1];
					zM = polygons->m[2][i+1];
				}
			}
			//c.green = (i * 50 + 50) % 255; //makes each surface visible
			//SHADING HERE-------------------------------------------------
			//ambient
			Ia.red = ambient.red * rcolor->r[0];
			Ia.green = ambient.green * rcolor->r[1];
			Ia.blue = ambient.blue * rcolor->r[2];

			//light vector
			light_v[0] = (xB + xM + xT) / 3.0 - point->l[0];
			light_v[1] = (yB + yM + yT) / 3.0 - point->l[1];
			light_v[2] = (zB + zM + zT) / 3.0 - point->l[2];

			//normalize
			normalize(light_v);
			normalize(normal);

			//diffuse
			Id.red = point->c[0] * rcolor->g[0] * calculate_dot(light_v, normal) * -1;
			Id.green = point->c[1] * rcolor->g[1] * calculate_dot(light_v, normal) * -1;
			Id. blue = point->c[2] * rcolor->g[2] * calculate_dot(light_v, normal) * -1;

			//specular
			reflect[0] = 2 * calculate_dot(light_v, normal) * normal[0] - light_v[0];
			reflect[1] = 2 * calculate_dot(light_v, normal) * normal[1] - light_v[1];
			reflect[2] = 2 * calculate_dot(light_v, normal) * normal[2] - light_v[2];
			Is.red = point->c[0] * rcolor->b[0] * calculate_dot(reflect, view) * calculate_dot(reflect, view) * calculate_dot(reflect, view);
			Is.green = point->c[1] * rcolor->b[1] * calculate_dot(reflect, view) * calculate_dot(reflect, view) * calculate_dot(reflect, view);
			Is.blue = point->c[2] * rcolor->b[2] * calculate_dot(reflect, view) * calculate_dot(reflect, view)  * calculate_dot(reflect, view);

			//add all I
			c.red = Ia.red + Id.red + Is.red;
			c.red = c.red>255?255:c.red;
			c.red = c.red<0?0:c.red;
			c.blue = Ia.blue + Id.blue + Is.blue;
			c.blue = c.blue>255?255:c.blue;
			c.blue = c.blue<0?0:c.blue;      
			c.green = Ia.green + Id.green + Is.green;
			c.green = c.green>255?255:c.green;
			c.green = c.green<0?0:c.green;
			//SHADING END---------------------------------------------------
			//3 outlines
			draw_line( polygons->m[0][i],
			 polygons->m[1][i],
			 polygons->m[2][i],
			 polygons->m[0][i+1],
			 polygons->m[1][i+1],
			 polygons->m[2][i+1],
			 s, c, zbuf);
			draw_line( polygons->m[0][i+1],
			 polygons->m[1][i+1],
			 polygons->m[2][i+1],
			 polygons->m[0][i+2],
			 polygons->m[1][i+2],
			 polygons->m[2][i+2],
			 s, c, zbuf);
			draw_line( polygons->m[0][i+2],
			 polygons->m[1][i+2],
			 polygons->m[2][i+2],
			 polygons->m[0][i],
			 polygons->m[1][i],
			 polygons->m[2][i],
			 s, c, zbuf);
			//fill in
			xL,xR = xB;
			zL,zR = zB;
			y = yB;
			//printf("draw_polygons\n");
			while(y<(int)yT){
			 if(y==(int)yB){
				 xL = xB;
				 zL = zB;
			 }
			 else{
				 d0 = (xT-xB)/(yT-yB);
				 xL = xB + d0*(y-yB);
				 d0 = (zT-zB)/(yT-yB);
				 zL = zB + d0*(y-yB);
		//printf("1%f %f\n",zL,zT);
			 }
			 if (y >= (int) yM){
				 if(y==(int)yM){
					 xR = xM;
					 zR = zM;
				 }
				 else{
					 d2 = (xT-xM)/(yT-yM);
					 xR = xM + d2*(y-yM);
					 d2 = (zT-zM)/(yT-yM);
					 zR = zM + d2*(y-yM);
			//printf("2%f %f\n",zR,zT);
				 }
			 }
			 else{
				 if(y==(int)yB){
					 xR = xB;
					 zR = zB;
				 }
				 else{
					 d1 = (xM-xB)/(yM-yB);
					 xR = xB + d1*(y-yB);
					 d1 = (zM-zB)/(yM-yB);
					 zR = zB + d1*(y-yB);
			//printf("3%f %f\n",xR,xM);
				 }
			 }
			 draw_line(xL,y,zL,xR,y,zR,s,c,zbuf);
			 y+=1;
		 }
		free(normal);
	 }
 }
}


/*======== void add_sphere() ==========
	Inputs:   struct matrix * points
						double cx
			double cy
			double r
			double step  
	Returns: 

	adds all the points for a sphere with center 
	(cx, cy) and radius r.

	should call generate_sphere to create the
	necessary points

	jdyrlandweaver
	====================*/

void add_sphere( struct matrix * points, 
		 double cx, double cy, double cz, double r, 
		 int step ) {

	struct matrix * temp;
	int lat, longt;
	int index;
	double x, y, z;
	int num_steps;
	
	num_steps = MAX_STEPS / step;

	temp = new_matrix( 4, num_steps * num_steps );
	//generate the points on the sphere
	generate_sphere( temp, cx, cy, cz, r, step );

	int latStop, longStop, latStart, longStart;
	latStart = 0;
	latStop = num_steps;
	longStart = 0;
	longStop = num_steps;
	
	for( lat = latStart; lat < latStop; lat++ ) {
		for ( longt = longStart; 2*longt < longStop; longt++ ) {
			
			index = lat * (num_steps) + longt;
			/*add_edge( points, temp->m[0][index],
		temp->m[1][index],
		temp->m[2][index],
		temp->m[0][index] + 1,
		temp->m[1][index] + 1,
		temp->m[2][index] );*/
			if( lat == latStop-1 ){
	int index2 = latStart * num_steps + longt;
	add_polygon( points,
				 temp->m[0][index], temp->m[1][index], temp->m[2][index],
				 temp->m[0][index2+1], temp->m[1][index2+1], temp->m[2][index2+1],
				 temp->m[0][index2], temp->m[1][index2], temp->m[2][index2]);
	add_polygon( points,
				 temp->m[0][index], temp->m[1][index], temp->m[2][index],
				 temp->m[0][index+1], temp->m[1][index+1], temp->m[2][index+1],
				 temp->m[0][index2+1], temp->m[1][index2+1], temp->m[2][index2+1]);
			}
			else{
	add_polygon( points, 
			 temp->m[0][index], temp->m[1][index], temp->m[2][index],
			 temp->m[0][index+num_steps+1], temp->m[1][index+num_steps+1], temp->m[2][index+num_steps+1],
			 temp->m[0][index+num_steps], temp->m[1][index+num_steps], temp->m[2][index+num_steps]);
	add_polygon( points,
			 temp->m[0][index], temp->m[1][index], temp->m[2][index],
			 temp->m[0][index+1], temp->m[1][index+1], temp->m[2][index+1],
			 temp->m[0][index+num_steps+1], temp->m[1][index+num_steps+1], temp->m[2][index+num_steps+1]);
			}
		}//end points only
		}
	free_matrix(temp);
}

/*======== void generate_sphere() ==========
	Inputs:   struct matrix * points
						double cx
			double cy
			double r
			double step  
	Returns: 

	Generates all the points along the surface of a 
	sphere with center (cx, cy) and radius r

	Adds these points to the matrix parameter

	03/22/12 11:30:26
	jdyrlandweaver
	====================*/

void generate_sphere( struct matrix * points, 
					double cx, double cy, double cz, double r, 
					int step ) {
	int circle, rotation;
	double x, y, z, circ, rot;

	int rotStart = step * 0;
	int rotStop = MAX_STEPS;
	int circStart = step * 0;
	int circStop = MAX_STEPS;
	
	for ( rotation = rotStart; rotation < rotStop; rotation += step ) {
		rot = (double)rotation / MAX_STEPS;
		for ( circle = circStart; circle < circStop; circle+= step ) {

			circ = (double)circle / MAX_STEPS;
			x = r * cos(2* M_PI * circ ) + cx;
			y = r * sin(2* M_PI * circ ) *
	cos( 2 * M_PI * rot ) + cy;
			z = r * sin(2* M_PI * circ ) *
	sin( 2 * M_PI * rot ) + cz;

			add_point( points, x, y, z);
		}
	}
}   
/*======== void add_torus() ==========
	Inputs:   struct matrix * points
						double cx
			double cy
			double r1
			double r2
			double step  
	Returns: 

	adds all the points required to make a torus
	with center (cx, cy) and radii r1 and r2.

	should call generate_torus to create the
	necessary points

	03/22/12 13:34:03
	jdyrlandweaver
	====================*/
void add_torus( struct matrix * points, 
		double cx, double cy, double cz, double r1, double r2, 
		int step ) {
	struct matrix * temp;
	int lat, longt;
	int index;
	int num_steps;
	
	num_steps = MAX_STEPS / step;

	temp = new_matrix( 4, num_steps * num_steps );
	//generate the points on the torus
	generate_torus( temp, cx, cy, cz, r1, r2, step );

	int latStop, longtStop, latStart, longStart;
	latStart = 0;
	longStart = 0;
	latStop = num_steps;
	longtStop = num_steps;
	for ( lat = latStart; lat < latStop; lat++ )
		for ( longt = longStart; longt < longtStop; longt++ ) {
			
			index = lat * num_steps + longt;
			
			if(lat == latStop-1){
				int index2 = latStart * num_steps + longt;
				if(longt == longtStop-1){
		add_polygon( points,
				 temp->m[0][index], temp->m[1][index], temp->m[2][index],
				 temp->m[0][index2-longt], temp->m[1][index2-longt], temp->m[2][index2-longt],
				 temp->m[0][index2], temp->m[1][index2], temp->m[2][index2]);
	add_polygon( points,
				 temp->m[0][index], temp->m[1][index], temp->m[2][index],
				 temp->m[0][index-num_steps+1], temp->m[1][index-num_steps+1], temp->m[2][index-num_steps+1],
				 temp->m[0][index2-longt], temp->m[1][index2-longt], temp->m[2][index2-longt]);
				}
				else{
				add_polygon( points,
							 temp->m[0][index], temp->m[1][index], temp->m[2][index],
							 temp->m[0][index2+1], temp->m[1][index2+1], temp->m[2][index2+1],
							 temp->m[0][index2], temp->m[1][index2], temp->m[2][index2]);
				add_polygon( points,
							 temp->m[0][index], temp->m[1][index], temp->m[2][index],
							 temp->m[0][index+1], temp->m[1][index+1], temp->m[2][index+1],
							 temp->m[0][index2+1], temp->m[1][index2+1], temp->m[2][index2+1]);
	}
			}
			else{
			if(longt == longtStop-1){
	add_polygon( points,
				 temp->m[0][index], temp->m[1][index], temp->m[2][index],
				 temp->m[0][index+1], temp->m[1][index+1], temp->m[2][index+1],
				 temp->m[0][index+num_steps], temp->m[1][index+num_steps], temp->m[2][index+num_steps]);
	add_polygon( points,
				 temp->m[0][index], temp->m[1][index], temp->m[2][index],
				 temp->m[0][index-longt], temp->m[1][index-longt], temp->m[2][index-longt],
				 temp->m[0][index+1], temp->m[1][index+1], temp->m[2][index+1]);
			}
			else{
	add_polygon( points, 
			 temp->m[0][index], temp->m[1][index], temp->m[2][index],
			 temp->m[0][index+num_steps+1], temp->m[1][index+num_steps+1], temp->m[2][index+num_steps+1],
			 temp->m[0][index+num_steps], temp->m[1][index+num_steps], temp->m[2][index+num_steps]);
	add_polygon( points,
			 temp->m[0][index], temp->m[1][index], temp->m[2][index],
			 temp->m[0][index+1], temp->m[1][index+1], temp->m[2][index+1],
			 temp->m[0][index+num_steps+1], temp->m[1][index+num_steps+1], temp->m[2][index+num_steps+1]);
			}
			}
		}//end points only
}

/*======== void generate_torus() ==========
	Inputs:   struct matrix * points
						double cx
			double cy
			double r
			double step  
	Returns: 

	Generates all the points along the surface of a 
	tarus with center (cx, cy) and radii r1 and r2

	Adds these points to the matrix parameter

	03/22/12 11:30:26
	jdyrlandweaver
	====================*/
void generate_torus( struct matrix * points, 
				 double cx, double cy, double cz, double r1, double r2, 
				 int step ) {

	double x, y, z, circ, rot;
	int circle, rotation;

	double rotStart = step * 0;
	double rotStop = MAX_STEPS;
	double circStart = step * 0;
	double circStop = MAX_STEPS;

	for ( rotation = rotStart; rotation < rotStop; rotation += step ) {

		rot = (double)rotation / MAX_STEPS;
		for ( circle = circStart; circle < circStop; circle+= step ) {

			circ = (double)circle / MAX_STEPS;
			x = cos( 2 * M_PI * rot ) *
	( r1 * cos( 2 * M_PI * circ ) + r2 ) + cx;
			y = r1 * sin( 2 * M_PI * circ ) + cy;
			z = sin( 2 * M_PI * rot ) *
	( r1 * cos( 2 * M_PI * circ ) + r2 ) + cz;

			add_point( points, x, y, z );
		}
	}
}

/*======== void add_box() ==========
	Inputs:   struct matrix * points
						double x
			double y
			double z
			double width
			double height
			double depth
	Returns: 

	add the points for a rectagular prism whose 
	upper-left corner is (x, y, z) with width, 
	height and depth dimensions.

	jdyrlandweaver
	====================*/
void add_box( struct matrix * polygons,
				double x, double y, double z,
				double width, double height, double depth ) {

	double x2, y2, z2;
	x2 = x + width;
	y2 = y - height;
	z2 = z - depth;
	//front
	add_polygon( polygons, 
				 x, y, z, 
				 x, y2, z,
				 x2, y2, z);
	add_polygon( polygons, 
				 x2, y2, z, 
				 x2, y, z,
				 x, y, z);
	//back
	add_polygon( polygons, 
				 x2, y, z2, 
				 x2, y2, z2,
				 x, y2, z2);
	add_polygon( polygons, 
				 x, y2, z2, 
				 x, y, z2,
				 x2, y, z2);
	//top
	add_polygon( polygons, 
				 x, y, z2, 
				 x, y, z,
				 x2, y, z);
	add_polygon( polygons, 
				 x2, y, z, 
				 x2, y, z2,
				 x, y, z2);
	//bottom
	add_polygon( polygons, 
				 x2, y2, z2, 
				 x2, y2, z,
				 x, y2, z);
	add_polygon( polygons, 
				 x, y2, z, 
				 x, y2, z2,
				 x2, y2, z2);
	//right side
	add_polygon( polygons, 
				 x2, y, z, 
				 x2, y2, z,
				 x2, y2, z2);
	add_polygon( polygons, 
				 x2, y2, z2, 
				 x2, y, z2,
				 x2, y, z);
	//left side
	add_polygon( polygons, 
				 x, y, z2, 
				 x, y2, z2,
				 x, y2, z);
	add_polygon( polygons, 
				 x, y2, z, 
				 x, y, z,
				 x, y, z2); 
}
	
/*======== void add_circle() ==========
	Inputs:   struct matrix * points
						double cx
			double cy
			double y
			double step  
	Returns: 


	03/16/12 19:53:52
	jdyrlandweaver
	====================*/
void add_circle( struct matrix * points, 
		 double cx, double cy, 
		 double r, double step ) {
	
	double x0, y0, x, y, t;
	
	x0 = cx + r;
	y0 = cy;

	for ( t = step; t <= 1; t+= step ) {
		
		x = r * cos( 2 * M_PI * t ) + cx;
		y = r * sin( 2 * M_PI * t ) + cy;
		
		add_edge( points, x0, y0, 0, x, y, 0 );
		x0 = x;
		y0 = y;
	}

	add_edge( points, x0, y0, 0, cx + r, cy, 0 );
}

/*======== void add_curve() ==========
Inputs:   struct matrix *points
				 double x0
				 double y0
				 double x1
				 double y1
				 double x2
				 double y2
				 double x3
				 double y3
				 double step
				 int type  
Returns: 

Adds the curve bounded by the 4 points passsed as parameters
of type specified in type (see matrix.h for curve type constants)
to the matrix points

03/16/12 15:24:25
jdyrlandweaver
====================*/
void add_curve( struct matrix *points, 
		double x0, double y0, 
		double x1, double y1, 
		double x2, double y2, 
		double x3, double y3, 
		double step, int type ) {

	double x, y, t;
	struct matrix * xcoefs;
	struct matrix * ycoefs;
	
	//generate the coeficients
	if ( type == BEZIER_MODE ) {
		ycoefs = generate_curve_coefs(y0, y1, y2, y3, BEZIER_MODE);
		xcoefs = generate_curve_coefs(x0, x1, x2, x3, BEZIER_MODE);
	}

	else {
		xcoefs = generate_curve_coefs(x0, x1, x2, x3, HERMITE_MODE);
		ycoefs = generate_curve_coefs(y0, y1, y2, y3, HERMITE_MODE);
	}

	for (t=step; t <= 1; t+= step) {
		
		x = xcoefs->m[0][0] * t * t * t + xcoefs->m[1][0] * t * t
			+ xcoefs->m[2][0] * t + xcoefs->m[3][0];

		y = ycoefs->m[0][0] * t * t * t + ycoefs->m[1][0] * t * t
			+ ycoefs->m[2][0] * t + ycoefs->m[3][0];

		add_edge(points, x0, y0, 0, x, y, 0);
		x0 = x;
		y0 = y;
	}

	free_matrix(xcoefs);
	free_matrix(ycoefs);
}

/*======== void add_point() ==========
Inputs:   struct matrix * points
				 int x
				 int y
				 int z 
Returns: 
adds point (x, y, z) to points and increment points.lastcol
if points is full, should call grow on points
====================*/
void add_point( struct matrix * points, double x, double y, double z) {
	
	if ( points->lastcol == points->cols )
		grow_matrix( points, points->lastcol + 100 );

	points->m[0][points->lastcol] = x;
	points->m[1][points->lastcol] = y;
	points->m[2][points->lastcol] = z;
	points->m[3][points->lastcol] = 1;

	points->lastcol++;
}

/*======== void add_edge() ==========
Inputs:   struct matrix * points
					int x0, int y0, int z0, int x1, int y1, int z1
Returns: 
add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
should use add_point
====================*/
void add_edge( struct matrix * points, 
				 double x0, double y0, double z0, 
				 double x1, double y1, double z1) {
	add_point( points, x0, y0, z0 );
	add_point( points, x1, y1, z1 );
}

/*======== void draw_lines() ==========
Inputs:   struct matrix * points
				 screen s
				 color c 
Returns: 
Go through points 2 at a time and call draw_line to add that line
to the screen
====================*/
void draw_lines( struct matrix * points, screen s, color c, struct matrix * zbuf) {

	int i;
 
	if ( points->lastcol < 2 ) {
		
		printf("Need at least 2 points to draw a line!\n");
		return;
	}

	for ( i = 0; i < points->lastcol - 1; i+=2 ) {

		draw_line( points->m[0][i], points->m[1][i], points->m[2][i],
				 points->m[0][i+1], points->m[1][i+1], points->m[2][i+1], s, c, zbuf);
		//FOR DEMONSTRATION PURPOSES ONLY
		//draw extra pixels so points can actually be seen    
		/*
		draw_line( points->m[0][i]+1, points->m[1][i], 
				 points->m[0][i+1]+1, points->m[1][i+1], s, c);
		draw_line( points->m[0][i], points->m[1][i]+1, 
				 points->m[0][i+1], points->m[1][i+1]+1, s, c);
		draw_line( points->m[0][i]-1, points->m[1][i], 
				 points->m[0][i+1]-1, points->m[1][i+1], s, c);
		draw_line( points->m[0][i], points->m[1][i]-1, 
				 points->m[0][i+1], points->m[1][i+1]-1, s, c);
		draw_line( points->m[0][i]+1, points->m[1][i]+1, 
				 points->m[0][i+1]+1, points->m[1][i+1]+1, s, c);
		draw_line( points->m[0][i]-1, points->m[1][i]+1, 
				 points->m[0][i+1]-1, points->m[1][i+1]+1, s, c);
		draw_line( points->m[0][i]-1, points->m[1][i]-1, 
				 points->m[0][i+1]-1, points->m[1][i+1]-1, s, c);
		draw_line( points->m[0][i]+1, points->m[1][i]-1, 
				 points->m[0][i+1]+1, points->m[1][i+1]-1, s, c);
		*/
	} 	       
}


void draw_line(int x0, int y0, double z0, int x1, int y1, double z1, screen s, color c, struct matrix* zbuf) {
 
	int x, y, d, dx, dy;
	double z, dz, dist;
	x = x0;
	y = y0;
	z = z0;
	//printf("%d %d %d %d\n", x0, x1, y0, y1);
	//dist = sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) );
	//swap points so we're always drawing left to right
	if ( x0 > x1 ) {
		//printf("flippy\n");
		x = x1;
		y = y1;
		z = z1;
		x1 = x0;
		y1 = y0;
		z1 = z0;
		x0 = x;
		y0 = y;
		z0 = z;
	}

	//need to know dx and dy for this version
	dx = (x1 - x) * 2;
	dy = (y1 - y) * 2;
	//printf("%d, %d, %d, %d, %d, %d, %d, %d\n",x,x0, x1,y, y0, y1, dx, dy);
	//printf("draw_line\n");
	//positive slope: Octants 1, 2 (5 and 6)
	if ( dy == 0 && dx == 0 ){
		//plot(s,c,x,y,z,zbuf);
	}
	else if ( dy > 0 ) {

		//slope < 1: Octant 1 (5)
		if ( dx > dy ) {
			d = dy - ( dx / 2 );
	
			while ( x <= x1 ) {
	plot(s, c, x, y, z, zbuf);
	if ( d < 0 ) {
		x = x + 1;
		d = d + dy;
	}
	else {
		x = x + 1;
		y = y + 1;
		d = d + dy - dx;
	}
	
	z = z0 + ((double)x-x0)/(x1-x0)*(z1-z0);
	//printf("1%f %f\n",z,z1);
			}
		}

		//slope > 1: Octant 2 (6)
		else {
			d = ( dy / 2 ) - dx;
			while ( y <= y1 ) {
	plot(s, c, x, y, z, zbuf);
	if ( d > 0 ) {
		y = y + 1;
		d = d - dx;
	}
	else {
		y = y + 1;
		x = x + 1;
		d = d + dy - dx;
	}
	z = z0 + ((double)y-y0)/(y1-y0)*(z1-z0);
	//printf("2%f %f\n",z,z1);
			}
		}
	}

	//negative slope: Octants 7, 8 (3 and 4)
	else { 

		//slope > -1: Octant 8 (4)
		if ( dx > abs(dy) ) {

			d = dy + ( dx / 2 );
	
			while ( x <= x1 ) {

	plot(s, c, x, y, z, zbuf);
	if ( d > 0 ) {
		x = x + 1;
		d = d + dy;
	}
	else {
		x = x + 1;
		y = y - 1;
		d = d + dy + dx;
	}
	z = z0 + ((double)x-x0)/(x1-x0)*(z1-z0);
	//printf("3%f %f\n",z,z1);
			}
		}

		//slope < -1: Octant 7 (3)
		else {

			d =  (dy / 2) + dx;

			while ( y >= y1 ) {
	
	plot(s, c, x, y, z, zbuf);
	if ( d < 0 ) {
		y = y - 1;
		d = d + dx;
	}
	else {
		y = y - 1;
		x = x + 1;
		d = d + dy + dx;
	}
	z = z0 + ((double)y-y0)/(y1-y0)*(z1-z0);
	//printf("4%f %f\n",z,z1);
			}
		}
	}
}

