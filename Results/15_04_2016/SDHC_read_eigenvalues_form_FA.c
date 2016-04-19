#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define USAGE "./load_env.x file_eigenvalue_1 file_eigenvalue_2 file_eigenvalue_3 file_output"
/*
  USAGE: ./load_env.x file_eigenvalue_1 file_eigenvalue_2 file_eigenvalue_3 file_output
  AUTHOR: Jaime Forero-Romero j.e.forero.romero@gmail.com
  MODIFIED BY: Sergio Hernandez Charpak sergiocharpak@gmail.com
  DESCRIPTION:
  Loads a 3Dimensional grid into a 3D array.
  NOTES:
  n_x, n_y, n_z is the grid size
  x_0, y_0, z_0 is the posiion of the cell for the grid(i,j,k) element in kpc/h
  dx,dy,z is the grid size in kpc/h
   
  grid[n] runs over x,y,z in space, following fortran array conventions
  for memory handling
   
  i.e. a 3D point grid(i,j,k) corresponds to  grid[i + N_X * (j + N_Y * k)] where N_X is the grid dimension size in x and N_Y is the grid dimension size in y (i,j,k start at 0).

*/

//First try - a little barbaric

//Definition of float size constant
#define FLOAT double

FLOAT * form_grid_eigen_value (FILE * in);
void getGeneralParameters(FILE * in, long long *n_total_gen, int * n_x_gen, int * n_y_gen, int * n_z_gen, long long * n_nodes_gen, float * dx_gen, float * dy_gen, float * dz_gen, float * x_0_gen, float * y_0_gen, float * z_0_gen);
void write_grid_FA(FILE *out, FLOAT*grid_FA, long long* n_nodes, long long* n_total_gen, int* n_x_gen, int* n_y_gen, int* n_z_gen,  float* dx_gen, float* dy_gen, float* dz_gen, float* x_0_gen, float* y_0_gen, float* z_0_gen);

int main(int argc, char **argv){
FILE *in;
FILE *out;
//Our 4 grids
FLOAT *grid_1;
FLOAT *grid_2;
FLOAT *grid_3;
FLOAT *grid_FA;

//These are general, they are valid for all 3 grids and we are using them in the formation of the fourth one

long long i;
long long n_total_gen;
int n_x_gen, n_y_gen, n_z_gen;
long long n_nodes_gen;
float dx_gen, dy_gen, dz_gen, x_0_gen, y_0_gen, z_0_gen;



if(argc!=5){
fprintf(stderr, "USAGE: %s\n", USAGE);
exit(1);
}

/*
//Calculation of the min and max in the grid to get the interval of values for the eigen value.
min_val = 1.0E10;
max_val = -1.0E10;
for(i=0;i<n_nodes;i++){
if(grid[i]<min_val){
min_val = grid[i];
}
if(grid[i]>max_val){
max_val = grid[i];
}
}
fprintf(stdout, "min val %g\n", min_val);
fprintf(stdout, "max val %g\n", max_val);
*/

//Now we begin the real thing.
//We get the general parameters from the first file
fprintf(stderr,"Getting the general parameters \n");
if(!(in=fopen(argv[1], "r"))){
fprintf(stderr, "Problem opening file %s\n", argv[1]);
exit(1);
}
getGeneralParameters(in, &n_total_gen, &n_x_gen, &n_y_gen, &n_z_gen, &n_nodes_gen, &dx_gen, &dy_gen, &dz_gen, &x_0_gen, &y_0_gen, &z_0_gen);

//Prints to verify
fprintf(stderr, "Nx Ny Nz : %d %d %d %d\n", n_x_gen, n_y_gen, n_z_gen, n_nodes_gen);
fprintf(stderr, "x_0 y_0 z_0 : %g %g %g\n", x_0_gen, y_0_gen, z_0_gen);
fprintf(stderr, "dx dy dz : %g %g %g\n", dx_gen, dy_gen, dz_gen);


//For the first EigenValue
fprintf(stderr,"Begins for the first eigenvalue \n");
if(!(in=fopen(argv[1], "r"))){
fprintf(stderr, "Problem opening file %s\n", argv[1]);
exit(1);
}

grid_1 = form_grid_eigen_value(in);

//For the second EigenValue
fprintf(stderr,"Begins for the second eigenvalue \n");
if(!(in=fopen(argv[2], "r"))){
fprintf(stderr, "Problem opening file %s\n", argv[2]);
exit(1);
}

grid_2 = form_grid_eigen_value(in);

//For the third Eigen Value
fprintf(stderr,"Begins for the third eigenvalue \n");
if(!(in=fopen(argv[3], "r"))){
fprintf(stderr, "Problem opening file %s\n", argv[3]);
exit(1);
}

grid_3 = form_grid_eigen_value(in);

fprintf(stderr,"We are done forming the three grids \n Now we proceed in the calculation of the FA \n");

//Allocates the grid
if(!(grid_FA=malloc(n_total_gen * sizeof(FLOAT)))){
fprintf(stderr, "problem with array allocation\n");
exit(1);
}

// We don't care about the coordinates. We just want the same position for the grids, here we have i
//FA = (1/sqrt(3))*(sqrt( ( (lambda1-lamda3)^2 + (lambda2 - lambda3)^2 + (lambda1 - lambda2)^2  )/( lambda1^2 + lambda2^2 + lambda3^2  )  ))
for(i=0;i<n_nodes_gen;i++){
grid_FA[i] = (1.0 / sqrt(3.0)) * sqrt( ( pow( (grid_1[i] - grid_3[i]), 2.0) + pow( (grid_2[i] - grid_3[i]), 2.0) + pow( (grid_1 - grid_2), 2.0)  )/( pow(grid_1[i], 2.0) + pow(grid_2[i], 2.0) + pow(grid_3[i], 2.0)  ) ); 
}

fprintf(stderr,"We are done forming the FA grid \n");

// We now proceed to write the grid_FA into a text file.
if(!(out = fopen(argv[4], "w"))){
fprintf(stderr, "Problem opening file %s\n", argv[4]);
exit(1);
}
 write_grid_FA(out, grid_FA, &n_nodes_gen, &n_total_gen, &n_x_gen, &n_y_gen, &n_z_gen, &dx_gen, &dy_gen, &dz_gen, &x_0_gen, &y_0_gen, &z_0_gen);

return 0;
}


/*
 * Gets the general parameters from one of the files to form the 4th grid
 */
void getGeneralParameters(FILE *in, long long* n_total_gen, int* n_x_gen, int* n_y_gen, int* n_z_gen, long long* n_nodes_gen, float* dx_gen, float* dy_gen, float* dz_gen, float* x_0_gen, float* y_0_gen, float* z_0_gen){
  
int dumb;
char line[30];
long long n_total;
int n_x, n_y, n_z;
long long n_nodes;
float dx, dy, dz, x_0, y_0, z_0;

fread(&dumb,sizeof(int),1,in);
fread(line,sizeof(char)*30,1,in);
fread(&dumb,sizeof(int),1,in);

fread(&dumb,sizeof(int),1,in);
fread(&n_x,sizeof(int),1,in);    
fread(&n_y,sizeof(int),1,in);    
fread(&n_z,sizeof(int),1,in);    
fread(&n_nodes,sizeof(long long),1,in);    
fread(&x_0,sizeof(float),1,in);    
fread(&y_0,sizeof(float),1,in);    
fread(&z_0,sizeof(float),1,in);    
fread(&dx,sizeof(float),1,in);    
fread(&dy,sizeof(float),1,in);    
fread(&dz,sizeof(float),1,in);    
fread(&dumb,sizeof(int),1,in);

n_total = n_x * n_y * n_z;

//Keeps the parameters
*n_total_gen = n_total;
*n_x_gen = n_x;
*n_y_gen = n_y;
*n_z_gen = n_z;
*n_nodes_gen = n_nodes;
 *dx_gen = dx;
 *dy_gen = dy;
 *dz_gen = dz;
 *x_0_gen = x_0;
 *y_0_gen = y_0;
 *z_0_gen = z_0;

 //Closes the file
 fclose(in);
}


/*
 * Forms the grid for the file incoming.
 */
FLOAT * form_grid_eigen_value (FILE * in){

  FLOAT *grid;
  int dumb;
  char line[30];
  long long i;
  long long n_total;
  int n_x, n_y, n_z;
  long long n_nodes;
  float dx, dy, dz, x_0, y_0, z_0;
  FLOAT max_val, min_val;

  fread(&dumb,sizeof(int),1,in);
  fread(line,sizeof(char)*30,1,in);
  fread(&dumb,sizeof(int),1,in);

  fread(&dumb,sizeof(int),1,in);
  fread(&n_x,sizeof(int),1,in);    
  fread(&n_y,sizeof(int),1,in);    
  fread(&n_z,sizeof(int),1,in);    
  fread(&n_nodes,sizeof(long long),1,in);    
  fread(&x_0,sizeof(float),1,in);    
  fread(&y_0,sizeof(float),1,in);    
  fread(&z_0,sizeof(float),1,in);    
  fread(&dx,sizeof(float),1,in);    
  fread(&dy,sizeof(float),1,in);    
  fread(&dz,sizeof(float),1,in);    
  fread(&dumb,sizeof(int),1,in);

  //Allocates the grid
  n_total = n_x * n_y * n_z;
  if(!(grid=malloc(n_total * sizeof(FLOAT)))){
    fprintf(stderr, "problem with array allocation\n");
    exit(1);
  }
  
  //Reads and fills the grid
  fread(&dumb,sizeof(int),1,in);
  fread(&(grid[0]),sizeof(FLOAT), n_nodes, in);
  fread(&dumb,sizeof(int),1,in);  

  fprintf(stderr, "Nx Ny Nz : %d %d %d %d\n", n_x, n_y, n_z, n_nodes);
  fprintf(stderr, "x_0 y_0 z_0 : %g %g %g\n", x_0, y_0, z_0);
  fprintf(stderr, "dx dy dz : %g %g %g\n", dx, dy, dz);
  fclose(in);

  //First three results of the grid
  fprintf(stdout, "%g %g %g  \n", grid[0], grid[1], grid[2]);
  
  return grid;

}


/*
 * Writes the FA grid into a text formatted file.
 */
void write_grid_FA(FILE *archivo, FLOAT*grid_FA, long long* n_nodes, long long* n_total_gen, int* n_x_gen, int* n_y_gen, int* n_z_gen, float* dx_gen, float* dy_gen, float* dz_gen, float* x_0_gen, float* y_0_gen, float* z_0_gen)
{

  //We have a big grid of n_x *n_y *n_z. We need to pass on some of the parameters.
  //the first line will be those parameters
  //The format will be: 
  //First Row: n_total_gen n_x_gen n_y_gen n_z_gen n_nodes_gen dx_gen dy_gen dz_gen x_0_gen y_0_gen z_0_gen
  //Division between rows, "\n" Division between pages: "\n \n" 

  /*
  float x0;
  float y0;
  x0=x[0];
  y0=y[0];
  char bufX[20];
  char bufY[20];
  char nmx= x0 -'0';
  char nmy= y0 -'0';
  int i;
  sprintf(bufX, "%f", x0);
  sprintf(bufY, "%f", y0);
  char n1[50], n3[50], n2[50];
   strcpy(n1,  "poblaciones_");
   strcpy(n2, "_");
   strcpy(n3, ".dat");

   strcat(n1, bufX);
   strcat(n1, n2);
   strcat(n1, bufY);
   strcat(n1, n3);
  */

  
   //First Line, the parameters
   fprintf(archivo, "%d \t %d \t %d \t %d \t %d \t %f \t %f \t %f \t %f \t %f \t %f \n", *n_nodes,  *n_total_gen, *n_x_gen, *n_y_gen,  *n_z_gen, *dx_gen, *dy_gen,  *dz_gen,  *x_0_gen,  *y_0_gen,  *z_0_gen);
   int i,j,k;
   //Finished the first line
  for(i=0;i<*n_x_gen;i++){
    for(j=0;j<*n_y_gen;j++){
      for(k=0;k<*n_z_gen;k++){
	fprintf(archivo, "%f", grid_FA[i + *n_x_gen * (j + *n_y_gen * k)]);
	fprintf(archivo, "\t");
      }
	  fprintf(archivo, "\n");
    }
    fprintf(archivo, "\n");
  }
  fclose(archivo);
}
