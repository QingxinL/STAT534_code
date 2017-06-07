//
//  main.cpp
//  test
//
//  Created by Yuxuan Cheng on 6/7/17.
//  Copyright © 2017 Yuxuan Cheng. All rights reserved.
//

#include <iostream>

int main(int argc, const char * argv[]) {
    
    // insert code here...
    int a=0;
    int b=5;
     //out;
    for (int i=0; i<30; i++)
    {
        
        double out = rand()/double(RAND_MAX);
        //int out = (rand() % (b-a+1))+ a;
        std::cout <<out<< "\n";
    }
    return 0;
}
