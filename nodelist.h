#ifndef NODELIST_HEADER
#define NODELIST_HEADER

#include <stdint.h>
#include "bmat.h"
#include "node.h"

#define DEFAULT_NODELIST_INCREASE 20

typedef struct nodelist_ 
{
	Node ** nodes;
	int * numChildren;
	int numNodes;
	int curNode;
	Node ** idxToNode;
	Node * root;
} NodeList;

void NodeList_init(NodeList * nl, int numNodes);
void NodeList_free(NodeList * nl);
Node * NodeList_get_Node(NodeList * nl);
void NodeList_resize(NodeList * nl, int numNewNodes);
void NodeList_add_size(NodeList * nl, int numNewNodes);
void NodeList_set_root(NodeList * nl, Node * root);
void NodeList_create_phylogeny(NodeList * nl, BMat * bmat, int * lj);
void NodeList_print_recursive_(Node * node, FILE * output, int indentLevel);
void NodeList_print(NodeList * nl, FILE * output);
void NodeList_get_num_children_recursive_(Node * node, int * numChildren);
void NodeList_get_num_children(NodeList * nl);
void NodeList_get_idxToNode_recursive_(Node * node, Node ** idxToNode);
void NodeList_get_idxToNode(NodeList * nl);

#endif
