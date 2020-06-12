/*extern "C" {} *///it will instruct the compiler to expect C linkage for your C functions, not C++ linkage.
//#include <thrust/find.h>
//#include <thrust/device_vector.h>
//#include <thrust/count.h>
//#include <thrust/copy.h>
//#include <thrust/execution_policy.h>
//#include <thrust/device_free.h>
#include <stdio.h>
//#include "all_structures.h"
#include "all_structure_cuda.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>


#include<vector>
#include <chrono> 


#define THREADS_PER_BLOCK 1024 //we can change it

using namespace std;
using namespace std::chrono;

__device__ int changeFlag1;
__device__ int changeFlag2;


__global__ void initialize(int nodes, int src, RT_Vertex* SSSP, double dummy_dist, int degree_frame)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int d = degree_frame;
	int gridstride = gridDim.x * blockDim.x;

	for (int index = idx; index < nodes; index += gridstride)
	{
		if (index == src) {
			SSSP[index].Root = src;
			SSSP[index].Parent = src;
			SSSP[index].Dist = 0.0;
			SSSP[index].EDGwt = 0;
		}
		else
		{
			SSSP[index].Root = src;
			int temp = index - 1;
			SSSP[index].Parent = (temp) / d;//we get the floor value as Root stores int
			/*//option1: all edge weight distinct and greater than (max edge_wt) * (n-1)
			SSSP[index].EDGwt = dummy_dist + index;
			int temp_dist = dummy_dist + index;
			while (temp > 0)
			{
				temp = temp / d;
				if (temp > 0)
				{
					temp_dist = temp_dist + dummy_dist + temp;
				}

			}
			SSSP[index].Dist = temp_dist;
			d_UpdatedDist[index] = temp_dist;*/

			//option 2: //taking all edges equal and greater than (max edge_wt) * (n-1)
			SSSP[index].EDGwt = dummy_dist;
			double temp_dist = dummy_dist;
			while (temp > 0)
			{
				temp = temp / d;
				if (temp > 0)
				{
					temp_dist = temp_dist + dummy_dist;
				}

			}
			SSSP[index].Dist = temp_dist;

		}

	}
}



__global__ void updateParent(int numS, xEdge_cuda* allChange_cuda, RT_Vertex* SSSP, int* changeFlag1, int* changeFlag2)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int gridstride = gridDim.x * blockDim.x;

	for (int index = idx; index < numS; index += gridstride)
	{

		//get the edge
		int node_1 = allChange_cuda[index].node1;
		int node_2 = allChange_cuda[index].node2;
		double edgeWeight = allChange_cuda[index].edge_wt;
		//reset it to 0
		//int Edgedone = 0;

		///method 1
		//if (SSSP[node_2].Dist > SSSP[node_1].Dist + edgeWeight)
		//{
		//	//Edgedone = 1;
		//	if (SSSP[node_2].Parent != node_1) {
		//		//Update Parent and EdgeWt
		//		SSSP[node_2].Parent = node_1;
		//		SSSP[node_2].EDGwt = edgeWeight;
		//		SSSP[node_2].Dist = SSSP[node_1].Dist + edgeWeight;
		//		*changeFlag1 = 1; //every time node dist is updated, the flag becomes 1
		//		//printf("*changeFlag1: %d", *changeFlag1);
		//		*changeFlag2 = 1;
		//	}
		//}

		///method 1A
		if (SSSP[node_2].Dist > SSSP[node_1].Dist + edgeWeight)
		{
			//Edgedone = 1;
			
				//Update Parent and EdgeWt
				SSSP[node_2].Parent = node_1;
				SSSP[node_2].EDGwt = edgeWeight;
				SSSP[node_2].Dist = SSSP[node_1].Dist + edgeWeight;
				*changeFlag1 = 1; //every time node dist is updated, the flag becomes 1
				//printf("*changeFlag1: %d", *changeFlag1);
				*changeFlag2 = 1;
			
		}

		//method 2
		//if (SSSP[SSSP[node_2].Parent].Dist + SSSP[node_2].EDGwt> SSSP[node_1].Dist + edgeWeight)
		//{
		//	if (SSSP[node_2].Parent != node_1) {
		//				//Update Parent and EdgeWt
		//				SSSP[node_2].Parent = node_1;
		//				SSSP[node_2].EDGwt = edgeWeight;
		//				//SSSP[node_2].Dist = SSSP[node_1].Dist + edgeWeight;
		//				*changeFlag1 = 1; //every time node dist is updated, the flag becomes 1
		//				//printf("*changeFlag1: %d", *changeFlag1);
		//				*changeFlag2 = 1;
		//			}
		//}
	}
}



__global__ void checkDist(int nodes, RT_Vertex* SSSP, int* changeFlag1, int* changeFlag2)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int gridstride = gridDim.x * blockDim.x;

	for (int index = idx; index < nodes; index += gridstride)
	{
		int parent = SSSP[index].Parent;
		if (SSSP[index].Dist != SSSP[parent].Dist + SSSP[index].EDGwt)
		{
			SSSP[index].Dist = SSSP[parent].Dist + SSSP[index].EDGwt;
			*changeFlag2 = 1;
			*changeFlag1 = 1;
		//printf("index: %d, dist: %f , parent's dist: %f, edgeweight: %f\n", index, SSSP[index].Dist, SSSP[parent].Dist, SSSP[index].EDGwt);
		}
	}
}


__global__ void helperFunction(int total_block, int thread_per_block, int nodes, int src, RT_Vertex* SSSP, double dummy_dist, int degree_frame, int totalChange, xEdge_cuda* allChange_cuda)
{
	if (threadIdx.x == 0)
	{
		initialize << <total_block, thread_per_block >> > (nodes, src, SSSP, dummy_dist, degree_frame);
		changeFlag1 = 1;
		changeFlag2 = 1;
		int itr = 0, itr2 = 0;

		///Method 1
		for (int i = 0; i < degree_frame; i++)
		{
			updateParent << < total_block, thread_per_block >> > (totalChange, allChange_cuda, SSSP, &changeFlag1, &changeFlag2);
			//cudaDeviceSynchronize();
			//printf("after itr %d changeFlag1: %d", itr, changeFlag1);
			itr++;
		}
		cudaDeviceSynchronize();
		changeFlag1 = 1;
		changeFlag2 = 1;
		
		while (changeFlag2 == 1)
		{
			changeFlag2 = 0;
			checkDist << <total_block, thread_per_block >> > (nodes, SSSP, &changeFlag1, &changeFlag2);
			cudaDeviceSynchronize();
			itr2++;
		}
		//cudaDeviceSynchronize();

		while (changeFlag1 == 1)
		{
			changeFlag1 = 0;
			updateParent << < total_block, thread_per_block >> > (totalChange, allChange_cuda, SSSP, &changeFlag1, &changeFlag2);
			cudaDeviceSynchronize();
			itr++;
			if (changeFlag1 == 1 /**&& itr % *degree_frame == 0**/)
			{
				changeFlag2 = 1;
				//printf("itr value is: %d \n", itr);
			}
			while (changeFlag2 == 1)
			{
				changeFlag2 = 0;
				checkDist << <total_block, thread_per_block >> > (nodes, SSSP, &changeFlag1, &changeFlag2);
				//cudaDeviceSynchronize();
				itr2++;

			}
			cudaDeviceSynchronize();
		}

		///Method 1A
		//cudaStream_t stream1;
		//cudaStreamCreateWithFlags(&stream1, cudaStreamNonBlocking);
		//cudaStream_t stream2;
		//cudaStreamCreateWithFlags(&stream2, cudaStreamNonBlocking);
		//changeFlag1 = 1;
		//changeFlag2 = 1;
		//while (changeFlag1 == 1 || changeFlag2 == 1)
		//{
		//	changeFlag1 = 0;
		//	updateParent << < total_block, thread_per_block , 0 , stream1 >> > (totalChange, allChange_cuda, SSSP, &changeFlag1, &changeFlag2);
		//	//cudaDeviceSynchronize();
		//	//printf("after itr %d changeFlag1: %d", itr, changeFlag1);
		//	itr++;
		//	changeFlag2 = 0;
		//	checkDist << <total_block, thread_per_block, 0, stream2 >> > (nodes, SSSP, &changeFlag1, &changeFlag2);
		//	//cudaDeviceSynchronize();
		//	itr2++;
		//	cudaDeviceSynchronize();
		//	
		//}
		//
		//cudaStreamDestroy(stream1);
		//cudaStreamDestroy(stream2);





		///Method 2
		//for (int i = 0; i < degree_frame; i++)
		//{
		//	updateParent << < total_block, thread_per_block >> > (totalChange, allChange_cuda, SSSP, &changeFlag1, &changeFlag2);
		//	//cudaDeviceSynchronize();
		//	//printf("after itr %d changeFlag1: %d", itr, changeFlag1);
		//	itr++;
		//}
		//cudaDeviceSynchronize();
		//changeFlag1 = 1;
		//changeFlag2 = 1;
		//while (changeFlag1 == 1)
		//{
		//	while (changeFlag1 == 1)
		//	{
		//		changeFlag1 = 0;
		//		updateParent << < total_block, thread_per_block >> > (totalChange, allChange_cuda, SSSP, &changeFlag1, &changeFlag2);
		//		//cudaDeviceSynchronize();
		//		//printf("after itr %d changeFlag1: %d", itr, changeFlag1);
		//		itr++;
		//	}
		//	while (changeFlag2 == 1)
		//	{
		//		changeFlag2 = 0;
		//		checkDist << <total_block, thread_per_block >> > (nodes, SSSP, &changeFlag1, &changeFlag2);
		//		//cudaDeviceSynchronize();
		//		itr2++;
		//		//printf("itr2 value: %d changeflag2: %d", itr2, changeFlag2);
		//	}
		//	cudaDeviceSynchronize();
		//}
		printf("itr %d itr2: %d", itr, itr2);

	}
}

void edge_update(int* totalChange, xEdge_cuda* allChange_cuda, RT_Vertex* SSSP, int total_block, int thread_per_block, int* nodes, int* degree_frame);

/*
1st arg: original graph file name
2nd arg: no. of nodes
3rd arg: no. of edges
4th arg: max edge weight
5th arg: degree of frame
6th: total_block
7th: thread_per_block
*/
int main(int argc, char* argv[]) {

	int nodes, edges, max_edgewt, degree_frame, total_block, thread_per_block;
	double dummy_dist;
	cudaError_t cudaStatus;
	nodes = atoi(argv[2]);
	edges = atoi(argv[3]);
	max_edgewt = atoi(argv[4]);
	printf("max edge weight: %d \n", max_edgewt);
	degree_frame = atoi(argv[5]);
	printf("frame degree: %d \n", degree_frame);
	dummy_dist = max_edgewt * (nodes - 1) + 1; //should be greater than (max edge_wt) * (n-1)
	total_block = atoi(argv[6]);
	thread_per_block = atoi(argv[7]);

	/*dummy_dist = 1677721401;*/
	printf("dummy dist: %f \n", dummy_dist);
	/*** Read original Graph ***/
	//int* colStartPtr_R;
	//cudaStatus = cudaMallocManaged((void**)&colStartPtr_R, (nodes + 1) * sizeof(int)); //we take nodes +1 to store the start ptr of the first row 
	//if (cudaStatus != cudaSuccess) {
	//	fprintf(stderr, "cudaMalloc failed!");
	//	/*goto Error;*/
	//}
	//int total_adjmatrix_size_R = edges * 2; //e.g.= (0 1 wt1), (1 0 wt1) both are same edge, but both will be there
	//Colwt2* cuda_adjlist_full_R;
	//cudaStatus = cudaMallocManaged(&cuda_adjlist_full_R, total_adjmatrix_size_R * sizeof(Colwt2));
	//if (cudaStatus != cudaSuccess) {
	//	fprintf(stderr, "cudaMalloc failed!");
	//	/*goto Error;*/
	//}

	xEdge_cuda* allChange_cuda; //store all edges as edges to be inserted in the frame
	int totalChange = edges; //we consider number of edges = no. of edges to be inserted
	cudaStatus = cudaMallocManaged(&allChange_cuda, totalChange * sizeof(xEdge_cuda));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed at changeEdge");
		/*goto Error;*/
	}

	//use below for direct path
	/*string file1 = "C:\\Users\\khand\\Desktop\\PhD\\CUDA test\\Test\\test 1\\fullGraph.txt";
	char* cstr1 = &file1[0];
	readin_graphU(&R, nodes, cstr1);*/

	//use below code if we use pass file name as argument
	//readin_graphU(&R, nodes, argv[1]);


	//use below code to pass the file name as relative path.
	//**keep the files in the same folder
	//string file1 = "./fullGraph.txt";
	//char* cstr1 = &file1[0];
	//readin_graphU4(colStartPtr_R, cuda_adjlist_full_R, cstr1, &nodes); //when local file used

	readin_graphU4(argv[1], &nodes, allChange_cuda); //when cmd line arg used

	cout << "Reading graph data successful" << endl;


	//Initializing  Rooted Tree Frame
	RT_Vertex* SSSP;
	cudaStatus = cudaMallocManaged(&SSSP, nodes * sizeof(RT_Vertex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed at SSSP structure");
		/*goto Error;*/
	}

	int src = 0; //the source from which the paths are computed

	//Time calculation
	auto startTime = high_resolution_clock::now();
	helperFunction << < 1, 1 >> > (total_block, thread_per_block, nodes, src, SSSP, dummy_dist, degree_frame, totalChange, allChange_cuda);
	cudaDeviceSynchronize();
	//Time calculation
	auto stopTime = high_resolution_clock::now();
	// Time calculation
	auto duration = duration_cast<microseconds>(stopTime - startTime);
	cout << "Total time taken for GPU operations: "
		<< duration.count() << " microseconds" << endl;

	//Test code start
	cout << "SSSP" << endl;
	/*for (int i = 0; i < nodes; i++)
	{
		cout << "*******" << endl;
		cout << "node" << i << endl << "dist" << SSSP[i].Dist << endl << "parent" << SSSP[i].Parent << endl;
	}*/
	for (int i = 0; i < 5; i++)
	{
		cout << "*******" << endl;
		cout << "node" << i << endl << "dist" << SSSP[i].Dist << endl << "parent" << SSSP[i].Parent << endl;
	}
	cout << "*******success*******" << endl;

	//Test code end

	cudaFree(allChange_cuda);
	cudaFree(SSSP);
	return 0;
}

