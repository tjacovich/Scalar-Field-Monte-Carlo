#include <stdlib.h>//Standard C library
#include <math.h>//C math library

/*Jacovich 6/16/2016*/
/*This generates the change in action due to a given field point, as well as the original action*/

void deltaction(double* phi, double temp, double a, double* dS_n, double* dS_o, int dim_lat, int cnt, double lambda)
  {
    int x = cnt%dim_lat;
    int t = cnt/dim_lat;

    double nei = 0;
    int xp = (x+1)%dim_lat; nei += phi[xp+dim_lat*t];
    xp = (x-1+dim_lat)%dim_lat; nei += phi[xp+dim_lat*t];
    int tp = (t+1)%dim_lat; nei += phi[x+dim_lat*tp];
    tp = (t-1+dim_lat)%dim_lat; nei += phi[x+dim_lat*tp];

    *dS_n = a*temp*temp/2-nei*temp+(lambda/4)*pow(temp,4);
    temp = phi[cnt];
    *dS_o = a*temp*temp/2-nei*temp+(lambda/4)*pow(temp,4);
 }
