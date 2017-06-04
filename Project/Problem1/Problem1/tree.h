/*
   FILE: TREE.H
*/

#ifndef _BIN_TREE
#define _BIN_TREE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <iomanip>
#include <time.h>

typedef struct myNode* LPNode;
typedef struct myNode Node;

struct myNode
{
   double key;

   LPNode Left; //left subtree
   LPNode Right; //right subtree
};

LPNode MakeNewNode(double key);
LPNode treeInsertOne(LPNode Root,LPNode newnode);
void treeInsert(LPNode& Root,double key);
void InorderTreeWalk(LPNode Root,double* sortedvector,int& nmax);
void DeleteTree(LPNode Root);
void printTree(LPNode Root,char* filename);
void ReduceTreeWalk(LPNode &Root, LPNode &UpperNode,
                    double* sortedvector,
                    int& nmax);
#endif
