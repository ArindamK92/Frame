#include <stdio.h>
#include <iostream>
//#include<list>
#include<vector> 
#include <fstream> 
#include <sstream>
using namespace std;

struct Edge {
	//Edge related parameters
	int node1;
	int node2;
	double edge_wt;
	void clear()
	{}
};


int readin_graphU4(char* myfile, vector<Edge>* allEdges)
{
	FILE* graph_file;
	char line[128];

	graph_file = fopen(myfile, "r");
	int l = -1, x=0;
	int prev_n1 = 0, prev_n2 = 0, prev_wt = 0;
	int n1, n2, wt;

	while (fgets(line, 128, graph_file) != NULL)
	{
		//Read line
		sscanf(line, "%d %d %d", &n1, &n2, &wt);
		if (n1 == n2)
		{
			continue;
			x++;
		}
		if (prev_n1 == n1 && prev_n2 == n2 && prev_wt > wt)
		{
			allEdges->at(l).edge_wt = wt;
			x++;
			continue;

		}

		l++;
		Edge edge;
		edge.node1 = n1;
		edge.node2 = n2;
		edge.edge_wt = wt;
		allEdges->push_back(edge);
		prev_n1 = n1;
		prev_n2 = n2;
		prev_wt = wt;

	}//end of while
	fclose(graph_file);

	return x;
}

//use to remove duplicate and self loop from all edges
int main2(int argc, char* argv[]) {

	vector<Edge> allEdges; //store all edges as edges to be inserted in the frame

	//cout << "start"<<endl;
	int x = readin_graphU4(argv[1], &allEdges); //when cmd line arg used
	//cout << "duplicate edges or self loop removed: " << x << endl;

	for (int i=0; i< allEdges.size();i++)
	{
		cout << allEdges.at(i).node1 << " " << allEdges.at(i).node2 << " " << allEdges.at(i).edge_wt << endl;
	}
	//cout << "end" << endl;
}