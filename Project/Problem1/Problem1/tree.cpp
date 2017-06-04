#ifndef _BIN_TREE
#include "tree.h"
#endif

//this function creates a new node
//with a given key
LPNode MakeNewNode(double key)
{
  LPNode mynode = new Node;
  mynode->Left = NULL;
  mynode->Right = NULL;
  mynode->key = key;

  return(mynode);
}

//inserts a new node in the tree
LPNode treeInsertOne(LPNode Root,LPNode newnode)
{
   if(NULL==Root)
   {
      return(newnode);
   }

   //insert smaller numbers in the left subtree
   if(newnode->key<Root->key)
   {
      Root->Left = treeInsertOne(Root->Left,newnode);
   }
   else //insert the other numbers in the right subtree
   {
      Root->Right = treeInsertOne(Root->Right,newnode);
   }
   return(Root);
}

//inserts a number in an existing binary tree
void treeInsert(LPNode& Root,double key)
{
   LPNode newnode = MakeNewNode(key);
   Root = treeInsertOne(Root,newnode);
   return;
}

//this function traverses the binary tree in inorder
//and harvests the keys in a vector "sortedvector"
//the keys will appear in their increasing order
//"nmax" gives the number of keys harvested so far
void InorderTreeWalk(LPNode Root,
                     double* sortedvector,
                     int& nmax)
{
   if(NULL!=Root)
   {
     //first traverse the left subtree
      InorderTreeWalk(Root->Left,
                      sortedvector,
                      nmax);

      //then the root
      sortedvector[nmax] = Root->key;
      nmax++;

      //and the right subtree
      InorderTreeWalk(Root->Right,
                      sortedvector,
                      nmax);
   }
   return;  
}

//deletes the entire tree
void DeleteTree(LPNode Root)
{
   if(NULL!=Root)
   {
      //deletes the left tree first
      DeleteTree(Root->Left);
      //then the right tree
      DeleteTree(Root->Right);
      //then the root
      delete Root;
   }
   return;
}

//auxiliary recursive function called by "printTree"
void printTreeRec(LPNode Root,FILE* out)
{
   if(NULL!=Root->Left)
   {
      printTreeRec(Root->Left,out);
      fprintf(out,"%.3lf -> %.3lf [label = \"LEFT\"];\n",Root->key,Root->Left->key);
   }
 
   if(NULL!=Root->Right)
   {
      fprintf(out,"%.3lf -> %.3lf [label = \"RIGHT\"];\n",Root->key,Root->Right->key);
      printTreeRec(Root->Right,out);
   }

   return;
}

//saves the edges of the tree in Graphviz format
void printTree(LPNode Root,char* filename)
{
   FILE* out = fopen(filename,"w");

   if(NULL==out)
   {
      fprintf(stderr,"Cannot open file [%s]\n",filename);
      return;
   }
   
   fprintf(out,"digraph G {\n");
   printTreeRec(Root,out);
   fprintf(out,"}\n");
   fclose(out);

   return;
}

