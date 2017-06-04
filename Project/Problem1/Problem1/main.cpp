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

int main()
{
    int i;
    int n = 8;
    
    char InputFile[] = "somenumbers.txt";
    char OutputFile[] = "mybinarytree.dot";
    
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
    
    //sort in increasing order by traversing the tree in inorder
    double* vectorIncreasing = allocvector(n);
    
    i = 0;
    InorderTreeWalk(mytree,vectorIncreasing,i);
    
    printf("Vector sorted in increasing order\n");
    printvector(vectorIncreasing,n);
    printf("\n\n");
    
    //print the tree in Graphviz format
    printTree(mytree,OutputFile);
    
    //delete the tree
    DeleteTree(mytree);
    
    //free the memory
    freevector(vector);
    freevector(vectorIncreasing);
    
    return(1);
}
