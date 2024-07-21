#include "Graph.h"
int main(int argc, char *argv[]) {
	puts("\n-----------------------------------------------------------------------------------------");
	Graph *graph = new Graph(argv[1], atoi(argv[2]));
// #ifdef _TEST_
	printf("#Filename=%s\n#K=%d\n",argv[1], atoi(argv[2]));
// #endif
	graph->read(); 
	graph->search(); 
	graph->write(); // delete graph;
	puts("-----------------------------------------------------------------------------------------\n");
}