#ifndef NODELIST_HEADER
#define NODELIST_HEADER

#include <stdint.h>
#include "bmat.h"
#include "node.h"

#define DEFAULT_NODELIST_INCREASE 20

typedef struct nodelist_ 
{
	Node ** nodes;
	int32_t * numChildren;
	int32_t numNodes;
	int32_t curNode;
	Node ** idxToNode;
	Node * root;
} NodeList;

void NodeList_init(NodeList * nl, int32_t numNodes);
void NodeList_free(NodeList * nl);
Node * NodeList_get_Node(NodeList * nl);
void NodeList_resize(NodeList * nl, int32_t numNewNodes);
void NodeList_add_size(NodeList * nl, int32_t numNewNodes);
void NodeList_set_root(NodeList * nl, Node * root);
void NodeList_create_phylogeny(NodeList * nl, BMat * bmat, int32_t * lj);
void NodeList_print_recursive_(Node * node, FILE * output, int32_t indentLevel);
void NodeList_print(NodeList * nl, FILE * output);
void NodeList_get_num_children_recursive_(Node * node, int32_t * numChildren);
void NodeList_get_num_children(NodeList * nl);
void NodeList_get_idxToNode_recursive_(Node * node, Node ** idxToNode);
void NodeList_get_idxToNode(NodeList * nl);

#endif
