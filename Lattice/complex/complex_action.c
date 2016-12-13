#include <stdlib.h>//Standard C library
#include <math.h>//C math library

/*Jacovich 6/16/2016*/
/*This generates the change in action due to a given field point, as well as the original action*/

void deltaction(double* phiR, double* phiI, double tempR, double tempI, double a,  double* dS_nR, double* dS_nI, double* dS_oR, double* dS_oI, int dim_lat, int cnt, double lambda, double h)
  {

    int x = cnt%dim_lat;
    int t = cnt/dim_lat;

    double neiR = 0;
    double neiI = 0;
    int xp = (x+1)%dim_lat; neiR += phiR[xp+dim_lat*t]; neiI += phiI[xp+dim_lat*t];
    xp = (x-1+dim_lat)%dim_lat; neiR += phiR[xp+dim_lat*t]; neiI += phiI[xp+dim_lat*t];
    int tp = (t+1)%dim_lat; neiR += phiR[x+dim_lat*tp]; neiI += phiI[x+dim_lat*tp];
    tp = (t-1+dim_lat)%dim_lat; neiR += phiR[x+dim_lat*tp]; neiI += phiI[x+dim_lat*tp];
 
    *dS_nR = a*(tempR*tempR/*+tempI*tempI*/)/2-tempR*neiR/*-tempI*neiI*/+(lambda/4)*(tempR*tempR/*+tempI*tempI*/)*(tempR*tempR/*+tempI*tempI*/)+h*(tempR/*+tempI*/);

    tempR = phiR[cnt];
    tempI = phiI[cnt];

    *dS_oR = a*(tempR*tempR/*+tempI*tempI*/)/2-tempR*neiR/*-tempI*neiI*/+(lambda/4)*(tempR*tempR/*+tempI*tempI*/)*(tempR*tempR/*+tempI*tempI*/)+h*(tempR/*+tempI*/);
 } 
