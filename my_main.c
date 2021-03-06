/*========== my_main.c ==========
  This is the only file you need to modify in order
  to get a working mdl project (for now).
  my_main.c will serve as the interpreter for mdl.
  When an mdl script goes through a lexer and parser,
  the resulting operations will be in the array op[].
  Your job is to go through each entry in op and perform
  the required action from the list below:
  frames: set num_frames (in misc_headers.h) for animation
  basename: set name (in misc_headers.h) for animation
  vary: manipluate knob values between two given frames
  over a specified interval
  set: set a knob to a given value
  setknobs: set all knobs to a given value
  push: push a new origin matrix onto the origin stack
  pop: remove the top matrix on the origin stack
  move/scale/rotate: create a transformation matrix
  based on the provided values, then
  multiply the current top of the
  origins stack by it.
  box/sphere/torus: create a solid object based on the
  provided values. Store that in a
  temporary matrix, multiply it by the
  current top of the origins stack, then
  call draw_polygons.
  line: create a line based on the provided values. Store
  that in a temporary matrix, multiply it by the
  current top of the origins stack, then call draw_lines.
  save: call save_extension with the provided filename
  display: view the image live
  jdyrlandweaver
  =========================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "parser.h"
#include "symtab.h"
#include "y.tab.h"
#include "misc_headers.h"
#include "matrix.h"
#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "stack.h"
/*======== void first_pass()) ==========
  Inputs:
  Returns:
  Checks the op array for any animation commands
  (frames, basename, vary)
  Should set num_frames and basename if the frames
  or basename commands are present
  If vary is found, but frames is not, the entire
  program should exit.
  If frames is found, but basename is not, set name
  to some default value, and print out a message
  with the name being used.
  jdyrlandweaver
  ====================*/
  void first_pass() {
  	int i;
  	char frame = 0;
  	char basename = 0;
  	for (i=0;i<lastop;i++) {
  		switch (op[i].opcode) {
  			case FRAMES:
  			if (frame) {
  				printf("Frames already defined.\n");
  			}
  			else {
  				num_frames = op[i].op.frames.num_frames;
  				frame = 1;
  			}
  			break;
  			case BASENAME:
  			strcpy(name, op[i].op.basename.p->name);
  			basename = 1;
  			break;
  			case VARY:
  			if (!frame) {
  				printf("Frames not defined before vary. Now die\n");
  				exit(0);
  			}
  			break;
  		}
  	}
  	if (frame && !basename) {
  		printf("Frames found but basename not defined. Default basename Anim given\n");
  		strcpy(name, "Anim");
  	}
  }
/*======== struct vary_node ** second_pass()) ==========
  Inputs:
  Returns: An array of vary_node linked lists
  In order to set the knobs for animation, we need to keep
  a seaprate value for each knob for each frame. We can do
  this by using an array of linked lists. Each array index
  will correspond to a frame (eg. knobs[0] would be the first
  frame, knobs[2] would be the 3rd frame and so on).
  Each index should contain a linked list of vary_nodes, each
  node contains a knob name, a value, and a pointer to the
  next node.
  Go through the opcode array, and when you find vary, go
  from knobs[0] to knobs[frames-1] and add (or modify) the
  vary_node corresponding to the given knob with the
  appropirate value.
  05/17/12 09:29:31
  jdyrlandweaver
  ====================*/
  struct vary_node ** second_pass() {
  	struct vary_node ** knobs = (struct vary_node **)malloc(sizeof(struct vary_node *) * num_frames);
  	int i,j;
  	double start_frame, end_frame, start_val, end_val;
  	for (i=0;i<num_frames;i++)
  		knobs[i] = 0;
  	for (i=0;i<lastop;i++) {
  		switch (op[i].opcode) {
  			case VARY:
  			for (j=0;j<num_frames;j++) {
  				struct vary_node * temp = (struct vary_node *)malloc(sizeof(struct vary_node));
  				strcpy(temp->name, op[i].op.vary.p->name);
  				start_frame = op[i].op.vary.start_frame;
  				end_frame = op[i].op.vary.end_frame;
  				start_val = op[i].op.vary.start_val;
  				end_val = op[i].op.vary.end_val;
  				if (j < start_frame)
  					temp->value = start_val;
  				else if (j > end_frame)
  					temp->value = end_val;
  				else
  					temp->value = (end_val - start_val)/(end_frame - start_frame)*(j - start_frame) + start_val;
  				temp->next = 0;
  				struct vary_node * here = knobs[j];
  				if (!here)
  					knobs[j] = temp;
  				else {
  					while (here->next && strcmp(here->name, temp->name))
  						here = here->next;
	  // printf("frame %d name %s value %d\n", j, here->name, here->value);
  					if (!strcmp(here->name, temp->name)) {
  						if (j >= start_frame)
  							here->value = temp->value;
  					}
  					else
  						here->next = temp;
  				}
  			}
  			break;
  		}
  	}
  	return knobs;
  }
/*======== void print_knobs() ==========
  Inputs:
  Returns:
  Goes through symtab and display all the knobs and their
  currnt values
  jdyrlandweaver
  ====================*/
  void print_knobs() {
  	int i;
  	printf( "ID\tNAME\t\tTYPE\t\tVALUE\n" );
  	for ( i=0; i < lastsym; i++ ) {
  		if ( symtab[i].type == SYM_VALUE ) {
  			printf( "%d\t%s\t\t", i, symtab[i].name );
  			printf( "SYM_VALUE\t");
  			printf( "%6.2f\n", symtab[i].s.value);
  		}
  	}
  }
/*======== void my_main() ==========
  Inputs:
  Returns:
  This is the main engine of the interpreter, it should
  handle most of the commadns in mdl.
  If frames is not present in the source (and therefore
  num_frames is 1, then process_knobs should be called.
  If frames is present, the enitre op array must be
  applied frames time. At the end of each frame iteration
  save the current screen to a file named the
  provided basename plus a numeric string such that the
  files will be listed in order, then clear the screen and
  reset any other data structures that need it.
  Important note: you cannot just name your files in
  regular sequence, like pic0, pic1, pic2, pic3... if that
  is done, then pic1, pic10, pic11... will come before pic2
  and so on. In order to keep things clear, add leading 0s
  to the numeric portion of the name. If you use sprintf,
  you can use "%0xd" for this purpose. It will add at most
  x 0s in front of a number, if needed, so if used correctly,
  and x = 4, you would get numbers like 0001, 0002, 0011,
  0487
  05/17/12 09:41:35
  jdyrlandweaver
  ====================*/
  void my_main( int polygons ) {
  	int i, f, j;
  	double step;
  	double xval, yval, zval, knob_value;
  	struct matrix *transform;
  	struct matrix *tmp;
  	struct stack *s;
  	screen t;
  	color g;
  	struct vary_node **knobs;
  	struct vary_node *vn;
  	char frame_name[128];
  	num_frames = 1;
  	step = 5;
  	g.red = 0;
  	g.green = 0;
  	g.blue = 0;
	//int shading;
  	struct matrix *zbuffer;

  	struct constants *rcolor;
  	rcolor = (struct constants *)malloc(sizeof(struct constants));
  	color ambient;
      	//int l_index = 0;
  	first_pass();

  	if (num_frames > 1)
  		knobs = second_pass();

  	for (f=0;f<num_frames;f++) {
	  struct light **point;
	  point = (struct light **)malloc(10*sizeof(struct light*));
	  int l_index = 0;
	  int shading;
  		printf("Frame %d - Thank You Mr. DW for Everything! (Spamming so you can't ignore this!!!)\n",f);
  		s = new_stack();
  		tmp = new_matrix(4, 1000);
  		transform = new_matrix(4,4);
	//set up zbuffer
  		zbuffer = new_matrix(XRES,YRES);
  		for (i = 0; i < zbuffer->rows; i++)
  			for (j = 0; j < zbuffer->cols; j++)
  				zbuffer->m[i][j]=-99999999;

  			clear_screen( t );
  			if (num_frames > 1) {
  				vn = knobs[f];
  				while (vn) {
  					set_value(lookup_symbol(vn->name), vn->value);
  					vn = vn->next;
  				}
	  //print_knobs();
  			}
  			for (i=0;i<lastop;i++) {
  				switch (op[i].opcode) {
				case SHADING:
				  if(!strcmp(op[i].op.shading.p->name,"FLAT"))
				    shading = FLAT;
				  else if(!strcmp(op[i].op.shading.p->name,"GOURAUD")){
				    shading = GOURAUD;
				    //printf("growled\n");
				  }
				  else if(!strcmp(op[i].op.shading.p->name,"PHONG"))
				    shading = PHONG;
				  else
				    shading = 0;
				  break;
  				case CONSTANTS:
				  // printf("CONSTANTS\n");
          rcolor = lookup_symbol(op[i].op.constants.p->name)->s.c;
  				break;
  				case AMBIENT:
				  // printf("AMBIENT\n");
  				ambient.red = op[i].op.ambient.c[0];
  				ambient.green = op[i].op.ambient.c[1];
  				ambient.blue = op[i].op.ambient.c[2];
  				break;
  				case LIGHT:
				  // printf("LIGHT\n");
          point[l_index] = lookup_symbol(op[i].op.light.p->name)->s.l;
	  l_index++;
          break;
  				case SPHERE:
				  // printf("START\n");
  					add_sphere( tmp,op[i].op.sphere.d[0], //cx
					op[i].op.sphere.d[1], //cy
					op[i].op.sphere.d[2], //cz
					op[i].op.sphere.r,
					step);
	//apply the current top origin
  					matrix_mult( s->data[ s->top ], tmp );
  					draw_polygons( tmp, t, g, zbuffer, rcolor, ambient, point, shading );
  					tmp->lastcol = 0;
					// printf("DONE\n");
  				break;
  				case TORUS:
					add_torus( tmp, op[i].op.torus.d[0], //cx
		  			op[i].op.torus.d[1], //cy
		   			op[i].op.torus.d[2], //cz
		   			op[i].op.torus.r0,
		   			op[i].op.torus.r1,
		   			step);
					matrix_mult( s->data[ s->top ], tmp );
					draw_polygons( tmp, t, g, zbuffer, rcolor, ambient, point, shading );
					tmp->lastcol = 0;
				break;
				case BOX:
					add_box( tmp, op[i].op.box.d0[0],
					op[i].op.box.d0[1],
					op[i].op.box.d0[2],
					op[i].op.box.d1[0],
					op[i].op.box.d1[1],
					op[i].op.box.d1[2]);
					matrix_mult( s->data[ s->top ], tmp );
					draw_polygons( tmp, t, g, zbuffer, rcolor, ambient, point, shading );
					tmp->lastcol = 0;
				break;
				/*case LINE:
					add_edge( tmp, op[i].op.line.p0[0],
					op[i].op.line.p0[1],
					op[i].op.line.p0[1],
					op[i].op.line.p1[0],
					op[i].op.line.p1[1],
					op[i].op.line.p1[1]);
					draw_lines( tmp, t, g, zbuffer, rcolor, ambient, shading );
					tmp->lastcol = 0;
				break;*/
				case MOVE:
				  //get the factors
				  // printf("MOVE\n");
				xval = op[i].op.move.d[0];
				yval = op[i].op.move.d[1];
				zval = op[i].op.move.d[2];
				knob_value = 1;
				if (op[i].op.move.p)
					knob_value = lookup_symbol(op[i].op.move.p->name)->s.value;
				transform = make_translate( xval * knob_value, yval * knob_value, zval * knob_value);
	//multiply by the existing origin
				matrix_mult( s->data[ s->top ], transform );
	//put the new matrix on the top
				copy_matrix( transform, s->data[ s->top ] );
				free_matrix( transform );
				break;
				case SCALE:
				xval = op[i].op.scale.d[0];
				yval = op[i].op.scale.d[1];
				zval = op[i].op.scale.d[2];
				knob_value = 1;
				if (op[i].op.scale.p)
					knob_value = lookup_symbol(op[i].op.scale.p->name)->s.value;
				transform = make_scale( xval * knob_value, yval * knob_value, zval * knob_value);
				matrix_mult( s->data[ s->top ], transform );
	//put the new matrix on the top
				copy_matrix( transform, s->data[ s->top ] );
				free_matrix( transform );
				break;
				case ROTATE:
				xval = op[i].op.rotate.degrees * ( M_PI / 180 );
				knob_value = 1;
				if (op[i].op.rotate.p)
					knob_value = lookup_symbol(op[i].op.rotate.p->name)->s.value;
	//get the axis
				if ( op[i].op.rotate.axis == 0 )
					transform = make_rotX( xval * knob_value);
				else if ( op[i].op.rotate.axis == 1 )
					transform = make_rotY( xval * knob_value);
				else if ( op[i].op.rotate.axis == 2 )
					transform = make_rotZ( xval * knob_value);
				matrix_mult( s->data[ s->top ], transform );
	//put the new matrix on the top
				copy_matrix( transform, s->data[ s->top ] );
				free_matrix( transform );
				break;
				case PUSH:
				push( s );
				break;
				case POP:
				pop( s );
				break;
				case SAVE:
				save_extension( t, op[i].op.save.p->name );
				break;
				case DISPLAY:
				display( t );
				break;
				}
			}
			//printf("SAVING\n");
			if(num_frames>1){
			sprintf(frame_name, "./anim/%s%03d.png", name, f);
			save_extension(t, frame_name);
			// printf("AFTERSAVE\n");
			}
			
			//free_stack( s );
			//free_matrix( tmp );
			//free_matrix( transform );
		// printf("VERY END\n");
	}
}
