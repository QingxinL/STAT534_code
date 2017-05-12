#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

int main( int argc, const char *argv[] ) {
    double x = 5.0 ;
    double y = gsl_sf_bessel_J0( x ) ;
    printf("J0(%g) = % .18e\n", x, y ) ;
    return 0 ;
}