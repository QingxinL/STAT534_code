//
//  main.cpp
//  test3
//
//  Created by Yuxuan Cheng on 4/11/17.
//  Copyright © 2017 Yuxuan Cheng. All rights reserved.
//

#include <iostream>
#include <math.h>
int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    printf ( "floor of 2.7 is %.1lf\n", floor (2.7) );
    double sigma = 2.0;
    int rd=0, cd=-5;
    double temp = 1/(2*3.14*sigma*sigma)*exp(-(rd*rd+cd*cd)/(2*sigma*sigma));
    printf ("%.7f\n",temp);
    std::cout<<sqrt(3)<<"    ";
   // double** image1;
    
    return 0;
}
