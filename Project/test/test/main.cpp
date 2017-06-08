//
//  main.cpp
//  test
//
//  Created by Yuxuan Cheng on 6/7/17.
//  Copyright © 2017 Yuxuan Cheng. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

double getPi(int n, int i)
{
    double Pi;
    Pi = (tgamma(n+1)/(tgamma(i+1)*tgamma(n-i+1)));
    return (Pi);
}

// return the random int belong [a, b]
int randInt(int a, int b)
{
    
    int out;
    out = (rand() % (b-a+1))+ a;
    return out;
}

// return the number [0,1]
double rand01()
{
    
    double out = rand()/double(RAND_MAX);
    return (out);
}

double MCMC(int n, int i, int iter)
{
    double sumI = 0;
    double p_is, p_it;
    int i_t;
    i_t = randInt(0, n);
    int i_t1 = 0;
    int i_s = 0;
    
    if (i_t==i) sumI++;  // for I_t0
    
    for (int j=0; j<iter; j++)
    {
        i_s = randInt(0, n);
        p_is = getPi(n, i_s);   // get the P_i_*
        p_it = getPi(n, i_t);   // get the P_i_t
        
        if (p_is>=p_it)
            i_t1 = i_s;
        else if (rand01() < (p_is/p_it))
            i_t1 = i_s;
        else
            i_t1 = i_t;
        
        i_t = i_t1;
        
        if (i_t==i) sumI++;
        
    }
    
    double p_est = double(sumI)/iter;
    return (p_est);
}

int main(int argc, const char * argv[]) {
    
    // insert code here...
    int a=0;
    int b=5;
    int n=10;
    
    std::cout<<getPi(10,0)<<std::endl;
     //out;
    //srand(unsigned(time(NULL)));
    /*
    for (int i=0; i<10; i++)
    {
        
        double out = randInt(a,b);
        //double out = rand01();
        //std::cout<<rand01()<<std::endl;
        std::cout <<out<< "\n";
    }
    */
    //std::cout<<tgamma(5)<<std::endl;
    for (int i=0; i<=n; i++)
    {
        std::cout<<MCMC(n, i, 25000)<<std::endl;
    }
    return 0;
}
