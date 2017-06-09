//
//  main.cpp
//  Problem1
//
//  Created by Yuxuan Cheng on 6/3/17.
//  Copyright © 2017 Yuxuan Cheng. All rights reserved.
//

#include <iostream>


/*
 This program sorts the elements of a given vector.
 */

#include "vectors.h"
#include "tree.h"

//open a file and return the count of numbers
int countNumbers(char* filename)
{
    int count=0;
    FILE* in = fopen(filename,"r");
    
    if(NULL==in)
    {
        printf("Cannot open file [%s]\n",filename);
    }
    
    double s;
    while(fscanf(in,"%lf",&s)==1)
    {
        count++;
    }
    
    fclose(in);
    
    return count;
}

int main()
{
    int i;
    //int n = 8;
    
    char InputFile[] = "somenumbers.txt";
    char OutputFile[] = "mybinarytree.dot";
//    char OutputFileTemp[] = "mybinarytreeTemp.dot";
    
    int n = countNumbers(InputFile);
    printf("n= %d\n", n);
    
    //initialise a vector of length n
    double* vector = allocvector(n);
    
    readvector(n,vector,InputFile);
    
    printf("Original vector\n");
    printvector(vector,n);
    printf("\n\n");
    
    //create a binary tree
    LPNode mytree = MakeNewNode(vector[0]);
    
    //now add all the other elements of the vector
    for(i=1;i<n;i++)
    {
        treeInsert(mytree,vector[i]);
    }
    
    //print the tree in Graphviz format
    printTree(mytree,OutputFile);
    
    // sort in increasing order by finding the smallest node in the tree
    int count = 0;
    double* vectorIncreasingNew = allocvector(n);
    for(i=0;i<n;i++)
    {
        ReduceTreeWalk(mytree, mytree,
                       vectorIncreasingNew,
                       count);
        
 //       printTree(mytree,OutputFileTemp);
        
    }
    
    printf("Vector sorted in increasing order\n");
    printvector(vectorIncreasingNew,n);
    printf("\n\n");
    
    
    //delete the tree
    DeleteTree(mytree);
    
    //free the memory
    freevector(vector);
    freevector(vectorIncreasingNew);
    
    return(1);
}
