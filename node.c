#include <stdlib.h>
#include <stdint.h>
#include "definitions.h"
#include "node.h"


void Node_init(Node * node)
{
	node->size = 0;
	node->root = 0;
	node->complete = 0;
	node->maxSize = INIT_NUM_CHILDREN;
	node->children = (Node **)malloc(sizeof(Node *) * (size_t)(node->maxSize));
	CHECKPOINTER(node->children);
	node->parent = NULL;
	return;
}

// adds numNewNodes to the size of the node's children
void Node_add_size(Node * node, int numNewNodes)
{
	node->children = (Node **)realloc((void *)node->children, sizeof(Node *) * (size_t)(node->maxSize + numNewNodes));
	CHECKPOINTER(node->children);
	node->maxSize += numNewNodes;
	return;
}

// resets the size of the node's children to newSize
void Node_resize(Node * node, int newSize)
{
	node->children = (Node **)realloc((void *)node->children, sizeof(Node *) * (size_t)(newSize));
	CHECKPOINTER(node->children);
	node->maxSize = newSize;
	return;
}

void Node_add_relationship(Node * parent, Node * child, int mutIdx)
{
	if(parent->size+1 > parent->maxSize)
		Node_resize(parent, DEFAULT_NODE_INCREASE);
	parent->children[parent->size] = child;
	parent->size++;
	child->parent = parent;
	child->mut = mutIdx;
	return;
}

/*
void Node_insert_node(Node * parent, Node * child, Node * insertNode, int mutIdx)
{
	int i;
	for(i = 0; i < parent->size; i++)
	{
		if(parent->children[i] == child)
		{
			parent->children[i] = insertNode;
			Node_add_relationship(insertNode, child);
			return;
		}
	}
	fprintf(stderr, "Warning: parent does not have child in Node_insert_node...nothing done\n");
	return;
}
*/

void Node_reset(Node * node)
{
	node->mut = -1;
	node->root = 0;
	node->parent = NULL;
	node->size = INIT_NUM_CHILDREN;
	Node_resize(node, INIT_NUM_CHILDREN);	// (resets maxSize)
	node->complete = 0;
	return;
}

void Node_free(Node * node)
{
	free(node->children);
	return;
}
