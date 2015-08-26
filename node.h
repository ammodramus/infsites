#ifndef NODE_HEADER
#define NODE_HEADER

#include <stdint.h>

#define INIT_NUM_CHILDREN 10
#define DEFAULT_NODE_INCREASE 10

typedef struct node_
{
	int mut;
	int root;
	int size;
	int maxSize;
	int complete;
	struct node_ ** children;
	struct node_ * parent;
} Node;

void Node_init(Node * node);
void Node_free(Node * node);
void Node_add_size(Node * node, int numNewNodes);
void Node_resize(Node * node, int numNewNodes);
void Node_add_relationship(Node * node, Node * child, int mutIdx);
void Node_insert_node(Node * node, Node * child, Node * insertNode);
void Node_reset(Node * node);

#endif
