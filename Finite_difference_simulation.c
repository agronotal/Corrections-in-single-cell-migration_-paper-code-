/*

  This code simulates the dynamics of the reaction diffusion system of equations


*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include "mt19937-64.c"

#define PI    3.14159265358979323846
#define NSNAP 1e3 // Number of snapshots

// Compute the sign of an argument
int sign(double x){
  return (x > 0) ? 1 : (-1);
}

double forceprex(double s, double prex, double pak){
  return  s-prex*pak +((prex*prex)/(1+(prex*prex)));
}

double forcepak(double a, double epsilon, double prex, double pak){
  return  -epsilon*(pak-a*prex);
}

// Choose next and previous elements with periodic boundary condition
int NEXT(int k, int N, int n) {
	return (k + n) % N;
}

int PREV(int k, int N, int n) {
	return ((k - n) % N + N) % N;
}

// Compute derivatives
// https://en.wikipedia.org/wiki/Finite_difference_coefficient
// https://en.wikipedia.org/wiki/Discrete_Laplace_operator

double Diff(double *Field, int i, double dx, int N){
	double FieldNext, FieldPrev;
	
	// Compute finite difference
	FieldNext = Field[NEXT(i,N,1)]*.8 - Field[NEXT(i,N,2)]*.2 + Field[NEXT(i,N,3)]*4./105 - Field[NEXT(i,N,4)]/280.;
	FieldPrev = Field[PREV(i,N,1)]*.8 - Field[PREV(i,N,2)]*.2 + Field[PREV(i,N,3)]*4./105 - Field[PREV(i,N,4)]/280.;

	return (FieldNext - FieldPrev) / dx;
}

double Diff2(double *Field, int i, double dx2, int N){
	double FieldNext, FieldPrev;

	// Compute finite difference
	FieldNext = Field[NEXT(i,N,1)]*1.6 - Field[NEXT(i,N,2)]*.2 + Field[NEXT(i,N,3)]*8./315 - Field[NEXT(i,N,4)]/560.;
	FieldPrev = Field[PREV(i,N,1)]*1.6 - Field[PREV(i,N,2)]*.2 + Field[PREV(i,N,3)]*8./315 - Field[PREV(i,N,4)]/560.;

	return (FieldNext + FieldPrev - Field[i]*205./72) / dx2;
}


/* ---------------------------
    DECLARATION OF FUNCTIONS
   ------------------------------*/
// Store data
void Record(FILE* outputdata, double *Density, double Time, int N);

// Implement the dynamics
void MoveDensity(int N, double *Noise_prex, double *prex,double *pak, double *fprex, double *fpak,  double *tder_prex,double *tder_pak, double s,double a,double epsilon, double diffprex,double diffpak,double dt, double dx2, double Time[0], double L);
/* ---------------------------
    END DECLARATION OF FUNCTIONS
   ------------------------------*/

/* ---------------------------
    MAIN
   ------------------------------*/

int main(int argc, char *argv[]){ 
  /* 
     argc is the number of argument of the command line, and argv is a list
     of string containing the arguments. argv[0] is the name of the
     executable.
  */
  chdir("data");
  char params[] = "number inputs incorrect";

  // Check that the number of inputs is correct
  if(argc!=13){
    printf("%s\n",params);
    exit(1);
  }
  
  
  // DEFINITION OF VARIABLES  
  int i, j, k; // Iterator
	int  N; // Number of lattice sites
  double L; // System size
  double s; // 
  double a; // 
  double epsilon;
  double diffprex;
  double diffpak;
  double *Noise_prex; // Array of noise for prex
  double *prex; // Array of prex field
  double *pak; // Array of pak field
  double *fprex; // Array of Prex reaction term values
  double *fpak; // Array of Pak reaction term values
  double *tder_prex; // Array of time derivative prex
  double *tder_pak; // Array of time derivative pak
  double dt; // Time step
  double dx, dx2; // Lattice constant
  double Time[0]; // Current time
  double totaltime; // Total duration of the run (after initialisation)
  double initial_time; // Time for initiation
  double interval_record; // Time between two recordings
  double nb_record;  // Number of time steps of the next recording
  double prex0; // initial prex value
  long counter_snap; // Counter for statistics
  long seed; // Seed of random number generator
  time_t clock; // Time to measure duration of the simulation
  

  char foldername[200];
  
  sprintf(foldername, "s_%s", argv[2]);
  mkdir(foldername); // creats folder with the relevant s
  chdir(foldername); // changes directory to that folder
  
  // READING THE INPUT
  i = 1;

  L               = strtod(argv[i], NULL); i++;
  s               = strtod(argv[i], NULL); i++;
  a               = strtod(argv[i], NULL); i++;
  epsilon         = strtod(argv[i], NULL); i++;
  diffprex        = strtod(argv[i], NULL); i++;
  diffpak         = strtod(argv[i], NULL); i++;
  dt              = strtod(argv[i], NULL); i++;
  totaltime       = strtod(argv[i], NULL); i++;
  initial_time    = strtod(argv[i], NULL); i++;
  interval_record = strtod(argv[i], NULL); i++;
  seed            = strtol(argv[i], NULL, 10); i++;
  dx              = strtod(argv[i], NULL); i++;
  dx2             = dx*dx; 
 


  FILE *output; // File where parameters are stored
  FILE *output_snap_prex; // File where snapshots of rho is stored
  FILE *output_snap_pak; // File where snapshots of m is stored
  char name[200]; // Name of output file containing the initial data
  char name_snap_prex[200]; // Name of output file containing snapshots of rho
  char name_snap_pak[200]; // Name of output file containing snapshots of m
  

  sprintf(name, "%.1f", a);
  sprintf(name_snap_prex, "%s_snap_prex", name);
  sprintf(name_snap_pak, "%s_snap_pak", name);
  

  output_snap_prex  = fopen(name_snap_prex, "w");
  output_snap_pak   = fopen(name_snap_pak, "w");
  output            = fopen(name, "w");


  // Print the parameters in output  
  fprintf(output, "%lg\n", L);
  fprintf(output, "%lg\n", s);
  fprintf(output, "%lg\n", a);
  fprintf(output, "%lg\n", epsilon);
  fprintf(output, "%lg\n", diffprex);
  fprintf(output, "%lg\n", diffpak);
  fprintf(output, "%lg\n", dt);
  fprintf(output, "%lg\n", totaltime);
  fprintf(output, "%lg\n", initial_time);
  fprintf(output, "%lg\n", interval_record);
  fprintf(output, "%lg\n", dx);
  fflush(output);

  // INITIALISATION OF VARIABLES
  init_genrand64(seed);

  // Number of sites
  N = (int) floor(L/dx);

  // Time counters
  nb_record = 0;


  // Start the clock time
  clock = time(NULL);

  // Initialize arrays
  Noise_prex        = (double*) malloc(sizeof(double)*N);
  prex              = (double*) malloc(sizeof(double)*N);
  pak               = (double*) malloc(sizeof(double)*N);
  fprex             = (double*) malloc(sizeof(double)*N);
  fpak              = (double*) malloc(sizeof(double)*N);
  tder_prex         = (double*) malloc(sizeof(double)*N);
  tder_pak          = (double*) malloc(sizeof(double)*N);

  // Initialaize Prex with small fluctuations around a homogeneous steady state

  prex0=sqrt((1 + s - a + sqrt(pow((1 + s - a),2) + 4*a*s))/(2*a));
  for(k=0;k<N;k++)	prex[k] = (prex0)*(1+0.1*cos(2*PI*k/N)); 
  for(k=0;k<N;k++)	pak[k]=a*prex0*(1+0.1*cos(2*PI*k/N));
  printf("Initial deposition succeeded\n");
  printf("s=%.3lg\n", s);
  printf("a=%.3lg\n", a);

  counter_snap = NSNAP;
  counter_snap = floor(totaltime/interval_record);
  Time[0]      = -initial_time;

  // Run dynamics until the time iterator reaches the final time
  while(Time[0]<totaltime){

    // Sample realizations of the noise 
    double noiseamp = s*(0.001); // Noise Amplitude
    for(k=0;k<N;k++) Noise_prex[k] = sqrt(2*noiseamp/dx)*gasdevMT(); 


    // Move the particles according to the dynamics

    MoveDensity(N, Noise_prex, prex, pak, fprex, fpak, tder_prex, tder_pak, s, a, epsilon, diffprex,diffpak,dt, dx2, Time,L);

		if(Time[0]>=nb_record){
	 		// Print simulation progress
			printf("%.3lg\n", 1e2*Time[0]/totaltime);
      
			// Record snapshot of system
			if(counter_snap>0){
				Record(output_snap_prex, prex, Time[0], N);
                Record(output_snap_pak, pak, Time[0], N);
				counter_snap--;
			}

			// Increase record counter
			nb_record += interval_record;

		}
  }

  // Return the duration of the simulation

  printf("Simulation duration: %ld seconds\n", (long) time(NULL)-clock);
  fprintf(output, "%ld", (long) time(NULL)-clock);
  
  return 0;
}

/* ---------------------------
    END MAIN
   ------------------------------*/

/* ---------------------------
    DYNAMICS
   ------------------------------*/

void MoveDensity(int N, double *Noise_prex, double *prex,double *pak, double *fprex, double *fpak,  double *tder_prex,double *tder_pak, double s,double a,double epsilon, double diffprex,double diffpak,double dt, double dx2, double Time[0], double L){
  int k; // Iterator
  double dprex; // Dummy variables
  double dpak; // Dummy variables
  double max_amp = 1e-1; // Maximum force amplitude
  double dt_g, sdt_g; // Time step

  // Reset reaction terms and time derivatives
  memset(fprex, 0, N*sizeof(double));
  memset(fpak, 0, N*sizeof(double));
  memset(tder_prex, 0, N*sizeof(double));
  memset(tder_pak, 0, N*sizeof(double));
	
	// Compute reaction term prex
  for(k=0;k<N;k++) fprex[k] = forceprex(s, prex[k], pak[k]);

// Compute reaction term pak
  for(k=0;k<N;k++) fpak[k] = forcepak( a,  epsilon, prex[k], pak[k]);


	// Compute time der. 
  for(k=0;k<N;k++){
	tder_prex[k] = fprex[k]+diffprex*Diff2(prex, k,  dx2,  N);
    tder_pak[k]    = fpak[k]  +diffpak*Diff2(pak, k,  dx2,  N);
	}


// Update maximum tder amplitude
  for(k=0;k<N;k++){
	if(fabs(tder_prex[k]) > max_amp) max_amp = fabs(tder_prex[k]);
    if(fabs(tder_pak[k]) > max_amp) max_amp = fabs(tder_pak[k]);
	}

        
 	// Adaptative time stepping
	if(max_amp*dt < 1e-1)	dt_g = dt;
	else									dt_g = 1e-1/max_amp;
	sdt_g    = sqrt(dt_g);
	Time[0] += dt_g;

  // Update density coordinate at all lattice sites
  for(k=0;k<N;k++){
	  // Dynamics
    dprex = dt_g*tder_prex[k]+sdt_g*Noise_prex[k] ;
    dpak   = dt_g*tder_pak[k] ;
    // New density after one iteration
    prex[k] += dprex;
    pak[k] += dpak;
  }
}

/* ---------------------------
    END DYNAMICS
   ------------------------------*/

/* ---------------------------
    PRINT
   ------------------------------*/

// Record particle coordinates
void Record(FILE* outputdata, double *Density, double Time, int N){

  int i;
  
  for(i=0;i<N;i++){
    fprintf(outputdata, "%lg\t%lg\n", Density[i], Time);
  }

  fflush(outputdata);
}

/* ---------------------------
    END PRINT
   ------------------------------*/
