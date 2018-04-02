#include <iostream>
#include <cmath>
#include "TF1.h"


double gres_ofy0_1(double x){return  sqrt(1.992*abs(1.992)/(x*x)+0.5539*0.5539*pow(x,-0.7915)+0.02733*0.02733);}
double gres_ofy0_2(double x){return  sqrt(3.301*abs(3.301)/(x*x)+0.6863*0.6863*pow(x,-0.8636)+0.02909*0.02909);}
double gres_ofy0_3(double x){return  sqrt( 4.42*abs( 4.42)/(x*x)+0.7891*0.7891*pow(x,-0.9105)+0.03046*0.03046);}
double gres_ofy0_4(double x){return  sqrt(5.566*abs(5.566)/(x*x)+0.8104*0.8104*pow(x,-0.9153)+0.03034*0.03034);}
double gres_ofy0_5(double x){return  sqrt(6.782*abs(6.782)/(x*x)+0.8801*0.8801*pow(x,-0.9417)+0.03078*0.03078);}

double wRho1=0.07858, wRho2=0.42870, wRho3=0.36185, wRho4=0.11092, wRho5=0.01995;

double gres_ofy0_all(double x){return  sqrt( wRho1*pow(gres_ofy0_1(x),2) +  wRho2*pow(gres_ofy0_2(x),2) + wRho3*pow(gres_ofy0_3(x),2) + wRho4* pow(gres_ofy0_4(x),2) + wRho5*pow(gres_ofy0_5(x),2) );}

void jetmetRes(){

TF1 *fres_ofy0_all = new TF1("fres_ofy0_all","gres_ofy0_all(x)",56,3132);

// use this function..



}
