#include <stdio.h>
#include <iostream>
//#include<list>
#include<vector> 
#include <fstream> 
#include <sstream>
using namespace std;

/******* Network Structures *********/


struct xEdge_cuda {
	//Edge related parameters
	int node1;
	int node2;
	double edge_wt;
	//End of Edge related parameters
	int inst;

	void clear()
	{}
};


// Data Structure for each vertex in the rooted tree
struct RT_Vertex
{
	int Root;    //root fo the tree
	int Parent; //mark the parent in the tree
	double EDGwt; //mark weight of the edge
	double Dist;  //Distance from root
};
//The Rooted tree is a vector of structure RT_Vertex;

////functions////
//Assumes the all nodes present
//Node starts from 0
//Total number of vertices=nodes and are consecutively arranged
//reads only the edges in the edge list, does not reverse them to make undirected


void readin_graphU4(char* myfile, int* nodes, xEdge_cuda* allChange_cuda)
{
	FILE* graph_file;
	char line[128];

	graph_file = fopen(myfile, "r");
	int l = 0;
	int prev_node = 0;

	while (fgets(line, 128, graph_file) != NULL)
	{
		int n1, n2, wt;
		//Read line
		sscanf(line, "%d %d %d", &n1, &n2, &wt);
		allChange_cuda[l].node1 = n1;
		allChange_cuda[l].node2 = n2;
		allChange_cuda[l].edge_wt = wt;
		allChange_cuda[l].inst = 1;
		prev_node = n1;
		l++;

	}//end of while
	fclose(graph_file);

	return;
}