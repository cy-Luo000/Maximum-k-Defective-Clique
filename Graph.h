#ifndef _GRAPH_H_
#define _GRAPH_H_
#include "Timer.h"
#include "LinearHeap.h"
class Graph {
private:
	std::string dir; //input graph directory
	int n; //number of nodes of the graph
	int m; //number of edges of the graph
	int K; //the value of k in k-plex

	int *pstart; //offset of neighbors of nodes
	int *pend; //used in search
	int *pend_buf;
	int *edges; //adjacent ids of edges
	int maxDeg;
	std::vector<int> KDC;

	int *s_degree;
	int *s_pstart;
	int *s_pend;
	int *s_edges;
	int *s_peel_sequence;
	int *s_core;
	char *s_vis;
	ListLinearHeap *s_heap;

	int *s_edgelist_pointer;
	int *s_tri_cnt;
	int *s_edge_list;
	int *s_active_edgelist;
	char *s_deleted;

public:
	Graph(const char *_dir, const int _K) ;
	~Graph() ;

	void read() ;
	void read_graph() ;

	void write() ;
	void verify_kplex() ;

	void search() ;

private:
	void extract_subgraph(int u, int *ids, int &ids_n, int *rid, std::vector<std::pair<int,int> > &vp, char *exists, int *pstart, int *pend, int *edges, char *deleted, int *edgelist_pointer) ;
	void extract_graph(int n, int m, int *deg, int *ids, int &ids_n, int *rid, std::vector<std::pair<int,int> > &vp, char *exists, int *pstart, int *pend, int *edges, char *deleted, int *edgelist_pointer) ;
	void induceSubgraph(int u, int *ids, int &ids_n, int *rid, std::vector<std::pair<int,int> > &vp, int *Q, int* degree, char *exists, int *pend, char *deleted, int *edgelist_pointer) ;

	int degen(int n, int *peel_sequence, int *core, int *pstart, int *edges, int *degree, char *vis, ListLinearHeap *heap, bool output) ;
	void shrink_graph(int &n, int &m, int *peel_sequence, int *core, int *out_mapping, int *in_mapping, int *rid, int *pstart, int *edges) ;
	void oriented_triangle_counting(int n, int m, int *peel_sequence, int *pstart, int *pend, int *edges, int *tri_cnt, int *adj) ;
	void reorganize_oriented_graph(int n, int *tri_cnt, int *edge_list, int *pstart, int *pend, int *pend2, int *edges, int *edgelist_pointer, int *buf) ;
	int peeling(int critical_vertex, ListLinearHeap *linear_heap, int *Qv, int &Qv_n, int d_threshold, int *Qe, int t_threshold, int *tri_cnt, int *active_edgelist, int &active_edgelist_n, int *edge_list, int *edgelist_pointer, char *deleted, int *degree, int *pstart, int *pend, int *edges, char *exists) ;
	char find(int u, int w, int &b, int e, char *deleted, int &idx, int *edgelist_pointer, int *edges) ;

	// functions for subgraph processing
	void load_graph_from_edgelist(int _n, const std::vector<std::pair<int,int> > &edge_list, int &n, int &m, int *degree, int *pstart, int *edges) ;
	int subgraph_heuri(ui *ids, ui &_n, std::vector<std::pair<int,int> > &edge_list, ui *rid, ui *Qv, ui *Qe, char *exists);
};
#endif
