#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "definitions.h"
#include "node.h"
#include "nodelist.h"

void NodeList_init(NodeList * nl, int32_t numNodes)
{
	int32_t i;
	nl->nodes = (Node **)malloc(sizeof(Node *) * (size_t)numNodes);
	CHECKPOINTER(nl->nodes);
	for(i = 0; i < numNodes; i++)
	{
		nl->nodes[i] = (Node *)malloc(sizeof(Node));
		CHECKPOINTER(nl->nodes[i]);
		Node_init(nl->nodes[i]);
	}
	nl->numChildren = (int32_t *)malloc(sizeof(int32_t) * (size_t)(numNodes));
	CHECKPOINTER(nl->numChildren);
	nl->idxToNode = (Node **)malloc(sizeof(Node *) * (size_t)(numNodes));
	CHECKPOINTER(nl->idxToNode);
	nl->numNodes = numNodes;
	nl->curNode = 0;
	nl->root = NULL;
	return;
}

void NodeList_free(NodeList * nl)
{
	int32_t i;
	for(i = 0; i < nl->numNodes; i++)
	{
		/* this assumes that each nl->nodes[i] has a ->children
		 * that has been allocated */
		Node_free(nl->nodes[i]);	
		free(nl->nodes[i]);
	}
	free(nl->nodes);
	free(nl->numChildren);
	free(nl->idxToNode);
	return;
}


Node * NodeList_get_Node(NodeList * nl)
{
	//if(nl->curNode+1 >= nl->numNodes)
		//NodeList_add_size(nl, DEFAULT_NODELIST_INCREASE);
	return nl->nodes[nl->curNode++];
}

void NodeList_resize(NodeList * nl, int32_t newNumNodes)
{
	nl->nodes = (Node **)realloc((void *)nl->nodes, sizeof(Node *) * (size_t)(newNumNodes));
	CHECKPOINTER(nl->nodes);
	nl->numNodes = newNumNodes;
	return;
}

void NodeList_add_size(NodeList * nl, int32_t numNewNodes)
{
	int32_t i, oldNumNodes = nl->numNodes;
	NodeList_resize(nl, numNewNodes + nl->numNodes);
	for(i = oldNumNodes; i < oldNumNodes+numNewNodes; i++)
		Node_init(nl->nodes[i]);
	return;
}

void NodeList_set_root(NodeList * nl, Node * root)
{
	nl->root = root;
	root->root = 1;
	root->mut = 0;
	return;
}

void NodeList_create_phylogeny(NodeList * nl, BMat * bmat, int32_t * lj)
{
	int32_t j, mutIdx = 1, ncols = bmat->ncols;
   	Node * root;
	Node ** colNodes = (Node **)malloc(sizeof(Node *) * (size_t)ncols);
	NodeList_init(nl, ncols+1);
	root = NodeList_get_Node(nl);
	root->size = 0;
	NodeList_set_root(nl, root);
	for(j = 0; j < ncols; j++)
	{
		colNodes[j] = NodeList_get_Node(nl);
	}
	for(j = 0; j < ncols; j++)
	{
		if(lj[j] > 0)
			Node_add_relationship(colNodes[lj[j]-1], colNodes[j], mutIdx++);
		else
			Node_add_relationship(root, colNodes[j], mutIdx++);
	}
	free(colNodes);
	return;
}

void NodeList_print_recursive_(Node * node, FILE * output, int32_t indentLevel)
{
	int32_t i;
	for(i = 0; i < indentLevel; i++)	
		fprintf(output, "    ");
	fprintf(output, "node = %p, mut = %i, root = %i, size = %i, maxSize = %i, complete = %i, parent = %p\n", node, node->mut, node->root, node->size, node->maxSize, node->complete, node->parent);
	for(i = 0; i < node->size; i++)
		NodeList_print_recursive_(node->children[i], output, indentLevel+1);
	return;
}

// assumes nl->root has been set
void NodeList_print(NodeList * nl, FILE * output)
{
	NodeList_print_recursive_(nl->root, output, 0);
	return;
}

void NodeList_get_num_children_recursive_(Node * node, int32_t * numChildren)
{
	int32_t i;
	numChildren[node->mut] = node->size;
	for(i = 0; i < node->size; i++)
		NodeList_get_num_children_recursive_(node->children[i], numChildren);
	return;
}

// assumes numChildren is a int32_t vector malloc'ed of length nl->numNodes
void NodeList_get_num_children(NodeList * nl)
{
	NodeList_get_num_children_recursive_(nl->root, nl->numChildren);
	return;
}

void NodeList_get_idxToNode_recursive_(Node * node, Node ** idxToNode)
{
	int32_t i;
	idxToNode[node->mut] = node;
	for(i = 0; i < node->size; i++)
		NodeList_get_idxToNode_recursive_(node->children[i], idxToNode);
	return;
}

void NodeList_get_idxToNode(NodeList * nl)
{
	NodeList_get_idxToNode_recursive_(nl->root, nl->idxToNode);
	return;
}
