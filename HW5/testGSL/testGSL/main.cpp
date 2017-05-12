//
//  main.cpp
//  testGSL
//
//  Created by Yuxuan Cheng on 5/11/17.
//  Copyright © 2017 Yuxuan Cheng. All rights reserved.
//

//#include <iostream>
//
//int main(int argc, const char * argv[]) {
//    // insert code here...
//    std::cout << "Hello, World!\n";
//    return 0;
//}

#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

int main( int argc, const char *argv[] ) {
    double x = 5.0 ;
    double y = gsl_sf_bessel_J0( x ) ;
    printf("J0(%g) = % .18e\n", x, y ) ;
    return 0 ;
}
