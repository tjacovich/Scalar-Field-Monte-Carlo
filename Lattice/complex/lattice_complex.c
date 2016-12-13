#include <stdlib.h>//Standard C Library
#include <math.h>//C math Library
#include <stdio.h>//C I/O library
#include <time.h>//Time functions library
#include "deltaction.h"//Action function header file
#include "mt64.h"//MT64 random number generator header file
#include "printlattice.h"//function to print the contents of the arrays
#include "printtwopoint.h"//Prints out an average of the correlator for a given field congfiguration
#include <complex.h>//Complex Numbers Library

/*Jacovich 8/1/16*/
/*All issues with action and number generation fixed.  Currently working with symmetry breaking
in 1+1 dimensions.  Implementing a complex scalar field as well.*/

/*Jacovich 6/29/16*/
/*Changed the random number generator to rand() in an attempt to fix issues with convergence.  
Found convergence with specific step and srand() seed for a given lattice size.  Not sure if
that is what should be expected.*/

/*Jacovich 6/16/16*/
/*This is the main program for the lattice field in 1+1 dimensions.  Future versions
will support 3+1 dimensions.  Action files can be swapped out by use of various headers*/

int main(int argc, char *argv[])
  { 
    /*checks to make sure the proper arguments are specified*/
    if(argc!=4)
     {
       printf("The number of arguments is invalid. \n");
       printf("usage: lattice_main dim_lat h file_name \n");
       return 1;
     }
     
    long arg1 = strtol(argv[1], &argv[1],10);//Converts command line arg into lattice spacing
    double arg2 = strtod(argv[2], &argv[2]);//Converts command line arg into symmetry breaking term 
    int dim_lat = arg1; /*Lattice Dimension*/
    double M = .25; /*The dimensionless mass parameter This parameter is related to the mass by M = MASS*LatticeSpacing*/
 // double* twopoint = calloc(dim_lat, sizeof(double));/*Holds Two Point function values*/
    double a = (4-M*M); /*defines the diagonal term of the kinetic action matrix*/
    double h = arg2; /*Symmetry breaking term*/
    int No = 100000; /*Number of field configurations to attempt*/
    int cnt = 1;/*Used for printing things out*/
    double step = 1.;/*determines the step-size of the perturbation*/
    double compstep = 0;/*step of the complex portion*/
    int acc = 0;/*number of acccepted proposals*/
    int tot = 0;/*total number of proposals*/
    double p = 2*4*atan(1)/10*0;//Momentum value
    double lambda = 0.005;/*coupling constant parameter*/
    
    /*seeding rand() in an attempt to further randomize*/
    init_genrand64(0); 
    
    double*  phiR = calloc(dim_lat*dim_lat,sizeof(double));/*main field array*/   
    //double*  phiI = calloc(dim_lat*dim_lat,sizeof(double));/*main field array*/
    double zero = 0;
    /*Defines the initial configuration of the field*/
    for(int l=0; l<dim_lat*dim_lat; l++)
      {
        phiR[l] = 2*step*(drand48() - 1);
       /* phiI[l] = 2*compstep*(drand48()-1);*/
      }
   
    FILE *outputfile = fopen(argv[3],"w+");//opens the output file for writing.

    /*This loop initiates a perturbation of a single field point on the lattice,
    and then implements a Metropolis style accept reject of the change in the field*/ 
   
    double tempI;/*imaginary component of temporary field*/
    double tempR;/*real component of temporary field*/
    double* S_I = calloc(dim_lat*dim_lat, sizeof(double));/*imaginary action*/  
    for(int i=0; i<No; i++)
      {
         for(int j=0; j<dim_lat*dim_lat; j++)
           {
             /*counts the number of iterations that occur*/ 
             tot++;
             
             /*Perturbs  the field with a pseudorandom real number between -1 and 1*/
             tempR = phiR[j] + step*( 2*drand48() - 1 ); 
            /* tempI = phiI[j] + compstep*( 2*drand48() - 1 );*/

             /*Generates the action due to the new field point*/  
             double S_nR, S_nI, S_oR, S_oI;
             deltaction(phiR, &zero, tempR, zero, a, &S_nR, &S_nI, &S_oR,  &S_oI, dim_lat, j, lambda, h);
             /*Metropolis Accept/Reject*/
             if(creal((exp(-S_nR)) / (exp(-S_oR))) > drand48())
               {
	               phiR[j] = tempR; //writes real part to array
	              // phiI[j] = tempI;//writes imaginary part to array
                 acc++; //counts the number of accepted values
               }
             
            } 
           

        /*Calculates the two point function for the most recent field conigfuration*/
        if(i>1000)
          {
            /*prints the momemntum space corellator*/
          //twopoint(phi, dim_lat, p, twopoint, outputfile);

            /*prints the average value of each field configuration*/
            double avgR = 0;
            double avgI = 0;
            double avgSR = 0;
            double avgSI = 0;

            /*calculates the complex value of the action*/
            for(int z = 0; z<dim_lat*dim_lat; z++)
              {
               // int x = z%dim_lat;
               // int t = z/dim_lat;
               // double nr = 0;
               // double nI = 0;
               // int xp = (x+1)%dim_lat; nr +=phiR[xp+t*dim_lat]; nI += phiI[xp+t*dim_lat];
               // int xm = (x-1+dim_lat)%dim_lat; nr +=phiR[xm+t*dim_lat]; nI += phiI[xm+t*dim_lat]; 
               // int tp = (t+1)%dim_lat; nr +=phiR[x+tp*dim_lat]; nI += phiI[x+tp*dim_lat];
               // int tm = (t-1+dim_lat)%dim_lat; nr +=phiR[x+tm*dim_lat]; nI += phiI[x+tm*dim_lat];
                S_I[z] =(+0/*phiI[z]*/);

              }
            
            for(int w=0;w<dim_lat*dim_lat;w++)
              {
                complex double b = (phiR[w] /*+ I*phiI[w]*/)*cexp(-I*S_I[w]);
                avgR += phiR[w]/(dim_lat*dim_lat);//creal(b)/(dim_lat*dim_lat);
                avgI += cimag(b)/(dim_lat*dim_lat);
                avgSR += creal(cexp(-I*S_I[w]))/(dim_lat*dim_lat);
                avgSI += cimag(cexp(-I*S_I[w]))/(dim_lat*dim_lat);
              }
            fprintf(outputfile, "%lf %lf %lf %lf \n", avgR, avgI, avgSR, avgSI);
          }
      }
    
   /*Prints the percentage of proposals accepted*/
   double bob = ((double) acc/tot);
   printf("%lf \n",bob);
								
   /*Free the allocated memory on the heap*/
   free(phiR);
  // free(phiI);
   free(S_I);
   fclose(outputfile);
   return 0;
 }
