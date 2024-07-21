#pragma once
#include "LinearHeap.h"
#include "Utility.h"
#include "Timer.h"
#include "ColorPacker.h"
using namespace std;
// #define _TEST_
#define DPBOUND
#define _SECOND_ORDER_REDUCTION_
// #define _DBUG_
// #define _CHECK_COIN_
#define _REMOVE_EDGES_
// #define _CLIQUEC_


int max_P_end=0;
int maxME=0; // the max num of missing edges in ans
long long tree_cnt=0, szSum=0, MuExSum=0;
int subNum=0;
int P_UBMin=10000;
int MaxMuExNum=0;
long long colMuExSum=0;
int maxColMuNum=0;
double denSum=0.0;
long long denNum=0;
double avgSubDense=0.0;
int maxSubSz=0;
long long colorBndPrune=0;
long long boundPrune=0;
int maxColSz=0;


class SubSearcher;
class ExactSearcher;
#ifdef _DBUG_
bool vecCmp(vector<int> &vec1, vector<int> &vec2){
	bool flag=true;
	if(vec1.size()!=vec2.size()){
		printf("unequal length: vec1: %d; vec2: %d\n",vec1.size(),vec2.size());
		return false;
	}

	int len=vec1.size();
	sort(vec1.begin(), vec1.end());
	sort(vec2.begin(),vec2.end());
	for (int i = 0; i < len; i++){
		if(vec1[i]!=vec2[i]) {
			printf("unequal val: vec1[%d]: %d; vec2[%d]: %d\n",i,vec1[i],i,vec2[i]);
			flag=false;
			break;
		}
	}
	if(!flag){
		printf("vec1: ");
		for(auto u: vec1) printf("%d, ",u);
		printf("\n");
		printf("vec2: ");
		for(auto u: vec2) printf("%d, ",u);
		printf("\n");
	}
	return flag;
}
void vecPush(vector<int>& vec,int *array,int start,int end){
	for (int i = start; i < end; i++)
		vec.push_back(array[i]);
	
}
#endif
class SubSearcher {
public:
	ui n; //number of nodes of the graph
	ept m; //number of edges of the graph
	//mapping
	// vector<ui> CNei;
	// vector<ui> CNeiRemap;
	ept *pstart; //offset of neighbors of nodes
	ept *pend; //used in search
	ui *edges; //adjacent ids of edges

	ui *tmp_vs;
	ui *tmp_color;

	ept *pstart_o; //oriented graph
	ui *edges_o; //oriented graph

	unsigned char *matrix; //adjacency matrix of an induced subgraph

	ui *head; // size of max_core, used for sorting vertices w.r.t. #colors
	ui *next; // size of max_core, used for sorting vertices w.r.t. #colors
	ui *id;

	ui maxCore;
	ui max_depth;
	long long branches;

	ui *degree;
	ui *rid;
	char *vis;

	ui *current_clique;
	std::vector<ui> vs_buf;
	std::vector<ui> color_buf;

	ui *mapping;
	ui mapping_n;
	ui matrix_len;

	std::vector<ui> MC;
	std::vector<ui> contractions;
	std::vector<std::pair<ui, ui> > changes;
	std::vector<ui> del, degree_one, degree_two, degree_three;

	SubSearcher() ;
	~SubSearcher() ;

	void read(const char *_file) ;
	void search() ;
	void write() ;
	void release() ;
	void load(ExactSearcher* exactSearcher) ;
	void loadCN(ExactSearcher* exactSearcher);

private:
	ui coloring_csr(const ui *vs, const ui vs_size, const ui original_size, ui *color, char *vis, const ui start_idx, const ui start_color) ;
	ui degen(ui *seq, ui *core, ui *color) ;
	ui degeneracy_maximal_clique_adj_list(const ui current_clique_size, ui *vs, const ui vs_size, char *vis, ui *degree, ListLinearHeap *heap) ;
	ui ego_degen(const ui *peel_sequence, const ui *core, const ui *color, ui *local_UBs, const ui UB) ;

	// build the oriented graph and store in pstart_o and edges_o, mapping of ids is stored in out_mapping
	// also reduces the sizes of peel_sequence, core, and color
	void shrink_graph(ui *&peel_sequence, ui *&core, ui *&color, ui *&out_mapping, ept *&pstart_o, ui *&edges_o) ;

	void search_oriented(const ui *peel_sequence, const ui *core, const ui *color, const ui *local_UBs) ;

	void recursive_search_clique_color_with_kernelization(const ui level, ui clique_size, ui vs_begin, ui vs_size, char kernel) ;
	void recursive_search_clique_color_without_kernelization(const ui level, ui clique_size, ui vs_begin, ui vs_size) ;

	void heuristic_max_clique_max_degree(ui processed_threshold) ;

	ui color_bound(const ui *vs, const ui vs_size, const ui *color, char *vis);

	ui color_bound(const ui *vs, const ui vs_size, const ui *mapping, const ui *color, char *vis);

	//coloring a graph that is represented by matrix
	//return the number of colors used
	ui coloring_matrix(const ui *vs, const ui vs_size, ui *color, char *vis, const ui start_idx, const ui start_color);

	//coloring a graph that is represented by matrix, aiming to minimize the number of vertices with color >= threshold
	//return the number of colors used
	ui coloring_matrix_advanced(const ui *vs, const ui vs_size, ui *color, const ui start_color, const ui threshold);

	// construct induced subgraph from pstart_o and edges_o
	void construct_induced_subgraph(const ui *vs, const ui vs_size, char *vis, ui *degree, const ept *pstart_o, const ui *edges_o);

	// construct matrix from pstart_o and edges_o for vertices in vs
	void construct_matrix(ui *vs, const ui vs_size, ui *mapping, char *vis, ui *rdegree, const ept *pstart_o, const ui *edges_o);

	ui degeneracy_maximal_clique_matrix(ui current_clique_size, ui *vs, const ui vs_size, ui *degree, char heuristic_gen, char print);

	void degree_one_two_reduction_with_folding_matrix(ui &current_clique_size, ui *vs, ui &vs_size, ui *rdegree, std::vector<ui> &contractions);

	void put_into_one_vector(ui &current_clique_size, ui &i, ui *vs, ui &vs_size, bool check, const ui rd, ui *rid);

	void put_into_one_vector_eq(ui &current_clique_size, ui &i, ui *vs, ui &vs_size, bool check, const ui rd, ui *rid);

	void degree_one_two_three_reduction_with_folding_matrix(ui &current_clique_size, ui *vs, ui &vs_size, ui *rdegree, ui *rid);

	void get_higher_neighbors(const ui u, ui &vs_size, std::vector<ui> &vs_buf, std::vector<ui> &color_buf, const ept *pstart_o, const ui *edges_o);

	//greedily enlarge max_clique by including u if feasible
	char greedy_extend(const ui u, const char *vis);

	void kcore_reduction(ui *vs, ui &vs_size, char *vis, ui *degree, const ui K, ui *queue);

	char kernelization_color(ui &clique_size, ui *vs, ui &vs_size, ui *color);

	void move_min_cardinality_color_to_front(ui *vs, const ui vs_size, ui *color);

	void obtain_degrees(const ui *vs, const ui vs_size, ui *degree);

	char reduce(ui *vs, ui &vs_size, ui *color, const ui threshold, ui clique_size);

	void remove_vertex_idx(ui *vs, ui idx, ui &vs_size, ui *rdegree);

	void search_triangle_matrix(const ui *vs, const ui vs_size, ui &clique_size);

	void search_triangle_matrix_color(const ui *vs, const ui vs_size, const ui *color, ui &clique_size);

	ui split_vs(ui *vs, const ui vs_size, ui *color, const ui threshold);

	void store_a_larger_clique(const ui clique_size, const char *info, char print);
};

class ExactSearcher{
public:
	int n;
	int m;
	int* pstart;
	int* edges;
    int k; //total missing edges
    int LB; // the lowerbound
    int UB; // the upperbound
	vector<int> &KDC;

    int* PC;
    int* PC_rid;// to record the positon in PC of each vertex v of G
    int P_end;// the P is from PC[0] to PC[P_end-1], used to enumerate defective set
    int C_end;// the C is from PC[P_end] to PC[C_end-1], represent R_{d}
	int CX_end;//the CX is from PC[C_end] to PC[CX_end-1], represent R_{0}, is the candidate set for clique

    int* neiInP;// to record every vertex's degree in P
    int* neiInG;// to record every vertex's degree in G
	int* neiInC; // to record every vertex's degree in C
	bool* matrix;
    int MEInP;// the missing edges in P 
    int MEInC;// the missing edges in C
	int MEInCCX; //the missing edges in C∪CX
    int MEInG;// the missing edges in G
	int MEInPC;// the missing edges in P∪C
	int *colB;// the bucket of the vertex of each color
	int CNColNum;// the color number of common neighbors

	long long treeIdx;
	
	int *colorLabel;// color label, begin from 0
	int *colorUse;
	int *colorSz;
	int *colorSzCX;
	int *colorIdx;
	int *colorUseMtx;
	int *weight;
	int *weightBucket;
	vector<vector<int>> nonNeiInPB;//used to bucket-sort the vertices by nonNeiInP 
	bool *colCX; 
	int *colVecC;
	int *colVecCX;
	int colNum;
	//define P_UB
	// int P_UB;
	//vertex removing
	int *remove_level;
	int *include_level;
	int *peel_level;
	int *restore_set;
	//the first pair: vertex or edge, the second pair: the move track, the third ele: if adjacent
	vector<tuple<pair<int,int>, pair<int,int>, int>> move_record;
	queue<int> Qv;
	queue<pair<int,int>> Qex;
#ifdef _SECOND_ORDER_REDUCTION_
	bool enable;
	bool usecsr;
	bool usemtx;
	bool useNew;
	bool delEdge=true;
	vector<vector<int>> CmNei;
	vector<vector<bool>> MuEx;
	queue<pair<int,int>> Qe;

#endif
#ifdef DPBOUND
	int *dp;
	int *t;
	colorPacker *colPacker;
	int *colExNum;
	int totExNum;
#endif
	int control=0;
#ifdef _DBUG_
	vector<vector<int>> P_rcd;
	vector<vector<int>> C_rcd;
	vector<vector<int>> CX_rcd;
	vector<int> MEInPRcd;
	vector<int> MEInCRcd;
	vector<int> MEInGRcd;
	vector<int> MEInCCXRcd;
	vector<int> MEInPCRcd;
#endif
	SubSearcher* subSearcher;


	bool& isAdj(int x, int y);
   ExactSearcher(int _k, int _n, vector<pair<int,int>> &_vp, vector<int> &_KDC, int _UB);
   ~ExactSearcher();
    
	bool isInP(int u);
	bool isInC(int u);
	bool isInCX(int u);
	bool isInX(int u);
	int checkEdges();
	int checkConinfidence(int level);
#ifdef _SECOND_ORDER_REDUCTION_
	void CalcuCmNei();

	bool edgePrune(int u, int v);
	int CalcuConflict(int *colVec, int start, int end, int level);
#endif
	void CtoP(int u, int level);
	bool CtoPwithPrune(int u, int level);
	void PtoC(int u, int level);
	void PtoCX(int u, int level);
	void CXtoC(int u, int level);
	void CtoCX(int u, int level);
	void CtoX(int u, int level);
	void XtoC(int u, int level);
	void CXtoX(int u, int level);
	void XtoCX(int u, int level);
	void swapID(int i, int j);
	void DSetBranch(int level, int P_UB);
	void search();
	
	bool judgeDsetFull(int P_UB);

	// int candiChoose();
	bool vertexCmp(int u, int v);

	int bound(int neis);
	int candiBucket();
#ifdef DPBOUND
	// int& DP(int c,int r){ return dp[c*n+r];}
	// int& T(int c, int r){return t[c*n+r];}
	

	void PrintColoring();
	int complexBound();
	bool bound();
	int PackingBound();
	int ColorBound();
	int CLUB();
#endif
	void store(int newLB);

	void maxClique4CN();
	int chooseVertex();
	void degreeSumPrune(int level);
	void colorBasedPrune(int level);

	bool includeV(int level, int PUB);
	void peelVertices(int level);
	bool collectRemovedEV(int level);
	bool removeVE(int level);


	void removeNonCN(int level);
	void recover(int level, int oldSz);
	void printPos(int u);
	bool existSatV();
	bool isDSet();
	void checkColUseMtx();
	

};

bool& ExactSearcher::isAdj(int x, int y){return matrix[x*n+y];}
ExactSearcher::ExactSearcher(int _k, int _n, vector<pair<int,int>> &_vp, vector<int> &_KDC, int _UB):KDC(_KDC){
		n=_n;
		m=_vp.size();
		k=_k;
		LB=_KDC.size();
		// LB=k+2;
		UB=_UB;

		pstart=new int[n+1];
		edges=new int[m*2];
        PC=new int[n];
        PC_rid=new int[n];
        neiInP=new int[n];
        neiInG=new int[n];
		neiInC = new int[n];
		matrix=new bool[n*n];

		colorLabel=new int[n];
		colorUse = new int[n];
		colorSz = new int[n];
		colorIdx = new int [n];
		colorUseMtx = new int[n*n];
		colB = new int[n*n];
		colCX = new bool[n];
		colVecC = new int[n*n];
		colVecCX = new int[n*n];
		colorSzCX = new int[n];
		CNColNum=0;
		colNum=0;
		fill(colVecC, colVecC+n*n, -1);
		fill(colVecCX, colVecCX+n*n, -1);
		memset(colorSzCX, 0, n*sizeof(int));
		weight = new int[n];
		//removing vertices
		remove_level = new int[n];
		include_level=new int[n];
		peel_level=new int[n];
		restore_set = new int[n];
		// P_UB=2*k+1;//the upper bound of P is initialize as 2*k+1
		nonNeiInPB.resize(k+1);
#ifdef DPBOUND
		dp =new int[n*n];
		t = new int[n*n];
		// colVec.resize(n);
		colPacker=new colorPacker(n);
		colExNum=new int[n];
		weightBucket=new int[n+k];
#endif
		P_end=0; C_end=n;treeIdx=0;
		CX_end=C_end; totExNum=0;
		memset(matrix,0,sizeof(bool)*n*n);
		memset(neiInP,0,sizeof(int)*n);
		memset(colExNum,0,n*sizeof(int));
		fill(colorLabel,colorLabel+n,-1);
		fill(colorIdx,colorIdx+n,-1);
		fill(colorUseMtx, colorUseMtx+n*n,-1);
		fill(colB, colB+n*n, -1);
		fill(remove_level, remove_level+n, n);
		fill(peel_level, peel_level+n, n);
		fill(include_level, include_level+n, n);
		fill(restore_set, restore_set+n, 1);
		memset(colorSz,0,n*sizeof(int));
		memset(colorUse,0,n*sizeof(int));
		memset(weight, 0, sizeof(int)*n);
		memset(weightBucket,0, (n+k)*sizeof(int));
		memset(colCX, false, n*sizeof(bool));
#ifdef _SECOND_ORDER_REDUCTION_
		enable=false;
		usecsr=true;
		usemtx=false;
		useNew=false;
		usecsr=false;
		usemtx=true;
		useNew=true;
		enable=true;
		// delEdge=false;
		CmNei.resize(n);
		MuEx.resize(n);
		for (int i = 0; i < n; i++){
			CmNei[i].resize(n,0);
			MuEx[i].resize(n, false);
		}	

#endif
#ifdef DPBOUND
		memset(dp, 0, n*n*sizeof(int));
		memset(t, 0, n*n*sizeof(int));
#endif
		for(auto pr : _vp){
			int u=pr.first, v=pr.second; isAdj(u,v)=isAdj(v,u)=true;
		}
		int idx=0;
		for(int i=0;i<n;++i){
			pstart[i]=idx;	
			for(int j=0;j<n;++j){
				if(isAdj(i,j)) edges[idx++]=j;
			}
			neiInG[i]=idx-pstart[i];
			neiInC[i]=neiInG[i];
			PC[i]=i;
			PC_rid[i]=i;
		}
		pstart[n]=idx;
		MEInP=0;
		MEInC=n*(n-1)/2-m;
		MEInPC=MEInC;
		MEInCCX = MEInC;
		MEInG=n*(n-1)/2-m;

		subSearcher=new SubSearcher;

		if(enable){
			for (int i = P_end; i < CX_end; i++){
				int u=PC[i];
				for (int j = i+1; j < CX_end; j++){
					int v=PC[j];
					for (int k = P_end; k < CX_end; k++){
						int w=PC[k];
						if(isAdj(u,w) && isAdj(v,w)) CmNei[u][v]++, CmNei[v][u]++;
					}
					
				}
			}
		}
#ifdef _DBUG_
		P_rcd.resize(n);
		C_rcd.resize(n);
		CX_rcd.resize(n);

		MEInPRcd.resize(n,0);
		MEInCRcd.resize(n,0);
		MEInGRcd.resize(n,0);
		MEInCCXRcd.resize(n,0);
		MEInPCRcd.resize(n,0);
		
#endif
}
ExactSearcher::~ExactSearcher(){
	delete[] PC;
	delete[] PC_rid;
	delete[] neiInC;
	delete[] neiInP;
	delete[] neiInG;
	delete[] matrix;
	delete[] pstart;
	delete[] edges;

	delete[] colorLabel;
	delete[] colorUse;
	delete[] colorSz;
	delete[] colorIdx;
	delete[] colorUseMtx;
	delete[] colB;
	delete[] colCX;
	delete[] colVecC;
	delete[] colVecCX;
	delete[] colorSzCX;
	delete[] weight;
	if(weightBucket!=NULL){
		delete[] weightBucket;
		weightBucket=NULL;
	}
	if(!nonNeiInPB.empty()) nonNeiInPB.clear();
	//vertex removing
	delete[] remove_level;
	delete[] peel_level;
	delete[] include_level;
	while(!Qv.empty()) Qv.pop();
	while(!Qex.empty()) Qex.pop();
#ifdef _SECOND_ORDER_REDUCTION_
	CmNei.clear();
	MuEx.clear();
	while (!Qe.empty()) Qe.pop();
	
#endif
#ifdef DPBOUND
	// colVec.clear();
	delete[] dp;
	delete[] t;
	if(colPacker!=NULL){
		delete colPacker;
		colPacker=NULL;
	}
	if(colExNum!=NULL){
		delete[] colExNum;
		colExNum=NULL;
	}
#endif

#ifdef _DBUG_
	P_rcd.clear();
	C_rcd.clear();
	CX_rcd.clear();
	MEInPRcd.clear();
	MEInCRcd.clear();
	MEInGRcd.clear();
	MEInCCXRcd.clear();
	MEInPCRcd.clear();

#endif
	delete subSearcher;
}
bool ExactSearcher::isInP(int u){return PC_rid[u] >= 0 && PC_rid[u] < P_end;}
bool ExactSearcher::isInC(int u){return PC_rid[u] >= P_end && PC_rid[u] < C_end;}
bool ExactSearcher::isInCX(int u){return PC_rid[u] >= C_end && PC_rid[u] < CX_end;}
bool ExactSearcher::isInX(int u){return PC_rid[u] >= CX_end && PC_rid[u] < n;}
int ExactSearcher::checkConinfidence(int level){
	int tMEInP=0, tMEInC=0, tMEInCCX=0, tMEInPC=0, tMEInG=0;
	for (int id = 0; id < CX_end; id++){
		int u=PC[id];
		int u_neiInP=0,u_neiInC=0,u_neiInG=0;
		int uu_neiInP=0,uu_neiInC=0,uu_neiInG=0;
		for (int i = pstart[u]; i < pstart[u+1]; i++){
			int vv=edges[i];
 			if(!isAdj(u,vv)) continue;
			if(isInP(vv)) uu_neiInP++, uu_neiInG++;
			else if(isInC(vv)) uu_neiInC++, uu_neiInG++;
			else if(isInCX(vv)) uu_neiInG++;
		}
		
		for (int j = 0; j < CX_end; j++){
			int v=PC[j];
 			if(isAdj(u,v)){
				if(isInP(v) || isInC(v) || isInCX(v)) u_neiInG++;
				if(isInP(v)) u_neiInP++;
				else if(isInC(v)) u_neiInC++;
			}
		}
		assert(u_neiInG==uu_neiInG);
		assert(u_neiInP==uu_neiInP);
		assert(u_neiInC==uu_neiInC);
		if(u_neiInP!=neiInP[u]){
			printf("level: %d, neiInP error: u_neiInP--%d, neiInP--%d\n",level, u_neiInP, neiInP[u]);
			return -1;
		}
 		if(u_neiInC!=neiInC[u]){
			printf("level: %d, neiInC error: u: %d, u_neiInC--%d, neiInC--%d\n",level, u,u_neiInC, neiInC[u]);
			return -1;
		}
		if(u_neiInG!=neiInG[u]){
			printf("level: %d, neiInG error: u_neiInG--%d, neiInG--%d\n",level, u_neiInG, neiInG[u]);
			return -1;
		}
		if(isInP(u)){
			tMEInP+=(P_end-1-neiInP[u]);
			tMEInPC+=(C_end-1-neiInP[u]-neiInC[u]);
			tMEInG+=(CX_end-1-neiInG[u]); 
		}else if(isInC(u)){
			tMEInC+=(C_end-P_end-1-neiInC[u]);
			tMEInPC+=(C_end-1-neiInP[u]-neiInC[u]);
			tMEInCCX+=(CX_end-P_end-1-(neiInG[u]-neiInP[u]));
			tMEInG+=(CX_end-1-neiInG[u]);
		}else if(isInCX(u)){
			tMEInCCX+=(CX_end-P_end-1-(neiInG[u]-neiInP[u]));
			tMEInG+=(CX_end-1-neiInG[u]);
		}
		
	}
	tMEInP/=2, tMEInC/=2, tMEInPC/=2, tMEInCCX/=2, tMEInG/=2;
	if(tMEInP!=MEInP){
		printf("level: %d, MEInP error: tMEInP--%d, MEInP--%d\n",level, tMEInP, MEInP);
		return -1;
	}
	if(tMEInC!=MEInC){
		printf("level: %d, MEInC error: tMEInC--%d, MEInC--%d\n",level, tMEInC, MEInC);
		return -1;
	}
	if(tMEInPC!=MEInPC){
		printf("level: %d, MEInPC error: tMEInPC--%d, MEInPC--%d\n",level, tMEInPC, MEInPC);
		return -1;
	}
	if(tMEInCCX!=MEInCCX){
		printf("level: %d, MEInCCX error: tMEInCCX--%d, MEInCCX--%d\n",level, tMEInCCX, MEInCCX);
		return -1;
	}
	if(tMEInG!=MEInG){
		printf("level: %d, MEInG error: tMEInG--%d, MEInG--%d\n",level, tMEInG, MEInG);
		return -1;
	}
	if(enable){
		vector<int> uneis;
		int cn=0;
		for (int i = 0; i < CX_end; i++){
			int u=PC[i];
			for (int j = i+1; j < CX_end; j++){
				int v=PC[j];
				cn=0;
				for (int k = P_end; k < CX_end; k++){
					int w=PC[k];
					if (isAdj(u,w)&&isAdj(v,w)) cn++;
				}
				if(cn!=CmNei[u][v] || cn!=CmNei[v][u]){
					printf("level: %d, cn error: cn--%d, cmNei[%d][%d]--%d\n",level, cn, u,v,CmNei[u][v]);
					return -1;
				}
				assert(cn==CmNei[u][v]);
				assert(cn==CmNei[v][u]);
			}
			
		}
		
	}
	return 1;
}
int ExactSearcher::checkEdges(){
	int edgeSum=0;
	for (int i = P_end; i < CX_end; i++){
		int u = PC[i];
		edgeSum+=(neiInG[u]-neiInP[u]);
	}
	return edgeSum;
}
void ExactSearcher::CalcuCmNei(){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++)
			CmNei[i][j]=0;
	}
	
	for (int i = 0; i < CX_end; i++){
		int u=PC[i];
		for (int j = i+1; j < CX_end; j++){
			int v=PC[j];
			for (int k = P_end; k < CX_end; k++){
				int w=PC[k];
				if(isAdj(u,w) && isAdj(v,w)) CmNei[u][v]++, CmNei[v][u]++;
			}
		}
	}
	return;
}

bool ExactSearcher::edgePrune(int u, int v){
	int cn=0;
	int e=(isAdj(u,v)==true)?0:1;
	if(!enable){
		for (int i = P_end; i < CX_end; i++)	{
			int w = PC[i];
			// if(neiInP[w]==P_end && isAdj(u,w) && isAdj(v,w)) cn++;
			if(isAdj(u,w) && isAdj(v,w)) cn++;
		}
		assert(cn==CmNei[u][v]);
	}else{
		cn=CmNei[u][v];
	}
	// assert(cn==CmNei[u][v]);
	int ub=P_end+2+cn+(k-MEInP-(P_end-neiInP[u]+P_end-neiInP[v])-e);
	// ub=n;
	int Ruv=k-(MEInP+P_end-neiInP[u]+P_end-neiInP[v])-e;
	if(isInP(u)&&!isInP(v)){
		ub=P_end+1+cn+(k-MEInP-(P_end-neiInP[v]));
		// ub=n;
		Ruv=k-(MEInP+P_end-neiInP[v]);
	}else if(!isInP(u)&&isInP(v)){
		ub=P_end+1+cn+(k-MEInP-(P_end-neiInP[u]));
		// ub=n;
		Ruv=k-(MEInP+P_end-neiInP[u]);
	}else if(isInP(u) && isInP(v)){
		ub=P_end+cn+(k-MEInP);
		Ruv=k-MEInP;
	}
	if (Ruv<0 || ub<=LB)
		return true;
	return false;
}

int ExactSearcher::CalcuConflict(int *colVec, int start, int end, int level){
	int conflictSum=0;
	for (int i = start; i < end; i++){
		int u = colVec[start];
		for (int j = i+1; j < end; j++){
			int v = colVec[j];
			if(MuEx[u][v]) conflictSum++;
		}
	}
	return conflictSum;
}
bool ExactSearcher::CtoPwithPrune(int u, int level){
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CtoPwithPrune begin\n", subNum, level);
		exit(0);
	}
#endif
	swapID(PC_rid[u], P_end++);
	int nonadj_inP=P_end-1-neiInP[u];
	int nonadj_inC=(C_end-P_end)-neiInC[u];
	int nonadj_inCCX = (CX_end-P_end)-(neiInG[u]-neiInP[u]);
	//renew MEC2P
	// MEC2P+=nonadj_inC;
	//renew MEInP
	MEInP+=nonadj_inP;
	//renew MEInC
	MEInC-=nonadj_inC;
	//renew the MEInCCX
	MEInCCX-=nonadj_inCCX;
	if(usemtx){
		vector<int> unei;
		for (int i = 0; i < CX_end; i++){
			int v=PC[i];
			if(isAdj(u,v)) unei.push_back(v);
		}
		for (int i = 0; i < unei.size(); i++){
			int v=unei[i];
			neiInP[v]++, neiInC[v]--;
			if(enable){
				for (int j = i+1; j < unei.size(); j++){
					int w=unei[j];
					CmNei[v][w]--;
					CmNei[w][v]--;
				}
				
			}
		}
		unei.clear();
	}
	//delete the exclusive vertices
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CtoPwithPrune end\n",subNum, level);
		exit(0);
	}
#endif
	return false;
}
void ExactSearcher::CtoP(int u, int level){
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CtoP begin\n", subNum, level);
		exit(0);
	}
#endif
#ifdef _DBUG_
	int edgeSum=checkEdges();
	if(edgeSum+2*MEInCCX!=(CX_end-P_end)*(CX_end-P_end-1)){
		printf("check edge fail: in CtoP\n");
		printf("P_end: %d, C_end: %d, CX_end: %d, MEInCCX: %d, edgesum_CCX: %d\n", P_end, C_end, CX_end, MEInCCX, edgeSum);
		exit(0);
	}
	assert(MEInCCX*2<=(CX_end-P_end)*(CX_end-P_end-1));
	// printf("enter CtoP 1\n");
	if(!isInC(u)) printPos(u);
	assert(isInC(u));
#endif
	// printf("enter CtoP 2\n");
	swapID(PC_rid[u], P_end++);
	int nonadj_inP=P_end-1-neiInP[u];
	int nonadj_inC=(C_end-P_end)-neiInC[u];
	int nonadj_inCCX = (CX_end-P_end)-(neiInG[u]-neiInP[u]);
	//renew MEInP
	MEInP+=nonadj_inP;
	//renew MEInC
	MEInC-=nonadj_inC;
	//renew the MEInCCX
	MEInCCX-=nonadj_inCCX;
	if (usemtx){
		vector<int> unei;
		for (int i = 0; i < CX_end; i++){
			int v=PC[i];
			if(isAdj(u,v)) unei.push_back(v);
		}
		for (int i = 0; i < unei.size(); i++){
			int v=unei[i];
			neiInP[v]++, neiInC[v]--;
			if(enable){
				for (int j = i+1; j < unei.size(); j++){
					int w=unei[j];
					CmNei[v][w]--;
					CmNei[w][v]--;
				}
				
			}
		}
		unei.clear();
	}
	
#ifdef _DBUG_
	assert(isInP(u));
	if (MEInCCX*2>(CX_end-P_end)*(CX_end-P_end-1)){
		printf("MEInCCX: %d, CX_end: %d, P_end: %d\n", MEInCCX, CX_end, P_end);
	}
	assert(MEInCCX*2<=(CX_end-P_end)*(CX_end-P_end-1));
#endif
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CtoP end\n", subNum, level);
		exit(0);
	}
#endif
}
void ExactSearcher::PtoC(int u, int level){
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in PtoC begin\n", subNum, level);
		exit(0);
	}
#endif
#ifdef _DBUG_
	if(!isInP(u)) printPos(u);
	assert(isInP(u));
	int edgeSum=checkEdges();
	if(edgeSum+2*MEInCCX!=(CX_end-P_end)*(CX_end-P_end-1)){
		printf("check edge fail: in PtoC\n");
		printf("subNum: %d, P_end: %d, C_end: %d, CX_end: %d, MEInCCX: %d, edgesum_CCX: %d\n", subNum, P_end,C_end, CX_end, MEInCCX, edgeSum);
		exit(0);
	}
#endif
	swapID(PC_rid[u], --P_end);
	int nonadj_inP=P_end-neiInP[u];
	int nonadj_inC=C_end-P_end-1-neiInC[u];
	int nonadj_inCCX=CX_end-P_end-1-(neiInG[u]-neiInP[u]);
	MEInP-=nonadj_inP;
	MEInC+=nonadj_inC;
	MEInCCX+=nonadj_inCCX;
	if (usemtx){
		vector<int> unei;
		for (int i = 0; i < CX_end; i++){
			int v=PC[i];
			if(isAdj(u,v)) unei.push_back(v);
		}
		for (int i = 0; i < unei.size(); i++){
			int v=unei[i];
			neiInP[v]--, neiInC[v]++;
			if(enable){
				for (int j = i+1; j < unei.size(); j++){
					int w=unei[j];
					CmNei[v][w]++;
					CmNei[w][v]++;
				}
			}
		}
		unei.clear();
	}
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in PtoC end\n", subNum, level);
		exit(0);
	}
#endif
}
void ExactSearcher::PtoCX(int u, int level){
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in PtoCX begin\n", subNum, level);
		exit(0);
	}
#endif
#ifdef _DBUG_
	if(!isInP(u)) printPos(u);
	int edgeSum=checkEdges();
	if(edgeSum+2*MEInCCX!=(CX_end-P_end)*(CX_end-P_end-1)){
		printf("check edge fail: in PtoCX\n");
		printf("P_end: %d, C_end: %d, CX_end: %d, MEInCCX: %d, edgesum_CCX: %d\n", P_end, C_end, CX_end, MEInCCX, edgeSum);
		exit(0);
	}
	assert(MEInCCX*2<=(CX_end-P_end)*(CX_end-P_end-1));
	assert(isInP(u));
#endif
	swapID(PC_rid[u], --P_end);
	swapID(P_end, --C_end);
	int nonadj_inP = P_end-neiInP[u];
	int nonadj_inPC = C_end-(neiInP[u]+neiInC[u]);
	int nonadj_inCCX = (CX_end - P_end - 1) - (neiInG[u] - neiInP[u]);
	MEInP -= nonadj_inP;
	MEInPC -= nonadj_inPC;
	MEInCCX += nonadj_inCCX;
	if (usemtx){
		vector<int> unei;
		for (int i = 0; i < CX_end; i++){
			int v=PC[i];
			if(isAdj(u,v)) unei.push_back(v);
		}
		for (int i = 0; i < unei.size(); i++){
			int v=unei[i];
			neiInP[v]--;
			if(enable){
				for (int j = i+1; j < unei.size(); j++){
					int w=unei[j];
					CmNei[v][w]++;
					CmNei[w][v]++;
				}
				
			}
		}
		unei.clear();
	}
	
#ifdef _DBUG_
	assert(isInCX(u));
	assert(MEInCCX*2<=(CX_end-P_end)*(CX_end-P_end-1));
#endif
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in PtoCX end\n", subNum, level);
		exit(0);
	}
#endif
}
void ExactSearcher::CXtoC(int u, int level){
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CXtoC begin\n", subNum, level);
		exit(0);
	}
#endif
#ifdef _DBUG_
	if(!isInCX(u)) {
		printPos(u);
		for (int i = 0; i < n; i++){
			printf("%d: %d, ",i, PC[i]);
		}
		printf("\n");
	}
	assert(isInCX(u));
	int edgeSum=checkEdges();
	if(edgeSum+2*MEInCCX!=(CX_end-P_end)*(CX_end-P_end-1)){
		printf("check edge fail: in CXtoC, level: %d, sub: %d\n",level, (CX_end-P_end)*(CX_end-P_end-1)-(edgeSum+2*MEInCCX));
		printf("subNum: %d, P_end: %d, C_end: %d, CX_end: %d, MEInCCX: %d, edgesum_CCX: %d\n", subNum,P_end, C_end, CX_end, MEInCCX, edgeSum);
		exit(0);
	}
	assert(isInCX(u));
	assert(MEInCCX*2<=(CX_end-P_end)*(CX_end-P_end-1));
#endif
	swapID(PC_rid[u], C_end++);
	int nonadj_InC = (C_end - P_end - 1) - neiInC[u];
	int nonadj_InPC = (C_end - 1) - (neiInP[u] + neiInC[u]);
	MEInC += nonadj_InC;
	MEInPC += nonadj_InPC;
	if (usemtx){
		vector<int> unei;
		for (int i = 0; i < CX_end; i++){
			int v=PC[i];
			if(isAdj(u,v)) unei.push_back(v);
		}
		for (int i = 0; i < unei.size(); i++){
			int v=unei[i];
			neiInC[v]++;
		}
		unei.clear();
	}
#ifdef _DBUG_
	assert(isInC(u));
	assert(MEInCCX*2<=(CX_end-P_end)*(CX_end-P_end-1));
#endif
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CXtoC end\n", subNum, level);
		exit(0);
	}
#endif
}
void ExactSearcher::CtoCX(int u, int level){
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CtoCX begin\n", subNum, level);
		exit(0);
	}
#endif
#ifdef _DBUG_
	assert(isInC(u));
#endif
	swapID(PC_rid[u], --C_end);
	int nonadj_InC = (C_end - P_end) - neiInC[u];
	int nonadj_InPC = C_end - (neiInP[u]+neiInC[u]);
	MEInC += nonadj_InC;
	MEInPC += nonadj_InPC;
	if (usemtx){
		vector<int> unei;
		for (int i = 0; i < CX_end; i++){
			int v=PC[i];
			if(isAdj(u,v)) unei.push_back(v);
		}
		for (int i = 0; i < unei.size(); i++){
			int v=unei[i];
			neiInC[v]--;
		}
		unei.clear();
	}
	
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CtoCX end\n", subNum, level);
		exit(0);
	}
#endif
}
void ExactSearcher::CtoX(int u, int level){
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CtoX begin\n", subNum, level);
		exit(0);
	}
#endif
#ifdef _DBUG_
	if(!isInC(u)) printPos(u);
	assert(isInC(u));
#endif
	swapID(PC_rid[u], --C_end);
	swapID(C_end, --CX_end);
	int nonadj_inC = (C_end-P_end)-neiInC[u];
	int nonadj_inPC = C_end - (neiInP[u]+neiInC[u]);
	int nonadj_inG = CX_end - neiInG[u];
	int nonadj_inCCX = (CX_end - P_end) - (neiInG[u] - neiInP[u]);
	MEInC -= nonadj_inC;
	MEInPC -= nonadj_inPC;
	MEInG -= nonadj_inG;
	MEInCCX -= nonadj_inCCX;
	if(usemtx){
		vector<int> unei;
		for (int i = 0; i < CX_end; i++){
			int v=PC[i];
			if(isAdj(u,v)) unei.push_back(v);
		}
		for (int i = 0; i < unei.size(); i++){
			int v=unei[i];
			neiInC[v]--, neiInG[v]--;
			if(enable){
				for (int j = i+1; j < unei.size(); j++){
					int w=unei[j];
					CmNei[v][w]--;
					CmNei[w][v]--;
				}
				
			}
		}
		unei.clear();
	}
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CtoX end\n", subNum, level);
		exit(0);
	}
#endif
}
void ExactSearcher::XtoC(int u,int level){
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in XtoC begin\n", subNum, level);
		exit(0);
	}
#endif
#ifdef _DBUG_
	if(!isInX(u)) printPos(u);
	assert(isInX(u));
#endif
	swapID(PC_rid[u], CX_end++);
	swapID(CX_end-1, C_end++);
	int nonadj_inC = (C_end-P_end-1)- neiInC[u];
	int nonadj_inPC = C_end - 1 - (neiInP[u] + neiInC[u]);
	int nonadj_inG = CX_end - 1 - neiInG[u];
	int nonadj_inCCX = (CX_end - P_end -1) - (neiInG[u] - neiInP[u]);
	MEInC += nonadj_inC;
	MEInPC += nonadj_inPC;
	MEInG += nonadj_inG;
	MEInCCX += nonadj_inCCX;
	if (usemtx){
		vector<int> unei;
		for (int i = 0; i < CX_end; i++){
			int v=PC[i];
			if(isAdj(u,v)) unei.push_back(v);
		}
		for (int i = 0; i < unei.size(); i++){
			int v=unei[i];
			neiInC[v]++, neiInG[v]++;
			if(enable){
				for (int j = i+1; j < unei.size(); j++){
					int w=unei[j];
					CmNei[v][w]++;
					CmNei[w][v]++;
				}
				
			}
		}
		unei.clear();
	}
	
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in XtoC end, u: %d\n", subNum, level,u);
		exit(0);
	}
#endif
}
void ExactSearcher::CXtoX(int u, int level){
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CXtoX begin\n", subNum, level);
		exit(0);
	}
#endif
#ifdef _DBUG_
	if(!isInCX(u)) printPos(u);
	assert(isInCX(u));
#endif
	swapID(PC_rid[u], --CX_end);
	int nonadj_inG = CX_end - neiInG[u];
	int nonadj_inCCX = (CX_end - P_end) - (neiInG[u] - neiInP[u]);
	MEInG -= nonadj_inG;
	MEInCCX -= nonadj_inCCX;
	if (usemtx){
		vector<int> unei;
		for (int i = 0; i < CX_end; i++){
			int v=PC[i];
			if(isAdj(u,v)) unei.push_back(v);
		}
		for (int i = 0; i < unei.size(); i++){
			int v=unei[i];
			neiInG[v]--;
			if(enable){
				for (int j = i+1; j < unei.size(); j++){
					int w=unei[j];
					CmNei[v][w]--;
					CmNei[w][v]--;
				}
				
			}
		}
		unei.clear();
	}
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in CXtoX end\n", subNum, level);
		exit(0);
	}
#endif
}
void ExactSearcher::XtoCX(int u, int level){
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in XtoCX begin\n", subNum, level);
		exit(0);
	}
#endif
#ifdef _DBUG_
	if(!isInX(u)) printPos(u);
	assert(isInX(u));
#endif
	swapID(PC_rid[u], CX_end++);
	int nonadj_inG = CX_end - 1 - neiInG[u];
	int nonadj_inCCX = (CX_end - P_end - 1)-(neiInG[u] - neiInP[u]);
	MEInG += nonadj_inG;
	MEInCCX += nonadj_inCCX;
	if (usemtx){
		vector<int> unei;
		for (int i = 0; i < CX_end; i++){
			int v=PC[i];
			if(isAdj(u,v)) unei.push_back(v);
		}
		for (int i = 0; i < unei.size(); i++){
			int v=unei[i];
			neiInG[v]++;
			if(enable){
				for (int j = i+1; j < unei.size(); j++){
					int w=unei[j];
					CmNei[v][w]++;
					CmNei[w][v]++;
				}
				
			}
		}
		unei.clear();
	}
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in XtoCX end, u: %d\n", subNum, level,u);
		exit(0);
	}
#endif
}
void ExactSearcher::DSetBranch(int level, int P_UB){
	// printf("enter dsetbranch\n");
	int PUB = min(P_UB, P_end+2*(k-MEInP));
	bool useInclude=true;
	bool usePeel=false;
	int oldSz=move_record.size();
	bool pruned=false;
	int oldSz2;
	int oldLB=LB, oldMEInP=MEInP;
#ifdef _DBUG_
	P_rcd[level].clear();
	vecPush(P_rcd[level],PC,0,P_end);
	C_rcd[level].clear();
	vecPush(C_rcd[level],PC,P_end,C_end);
	CX_rcd[level].clear();
	vecPush(CX_rcd[level],PC,C_end,CX_end);
	MEInPRcd[level]=MEInP;
	MEInCRcd[level]=MEInC;
	MEInGRcd[level]=MEInG;
	MEInPCRcd[level]=MEInPC;
	MEInCCXRcd[level]=MEInCCX;
#endif

	P_UBMin=min(P_UBMin, PUB);
	max_P_end = max(max_P_end,P_end);
	treeIdx++;
	if(useInclude){
		if(includeV(level,PUB)) {    //CtoP
			if(!isDSet()) goto REC;
			removeNonCN(level); //CtoX, CXtoX
			maxClique4CN();
			goto REC;}
	}
	PUB=min(PUB, P_end+2*(k-MEInP));
	// CalcuCmNei();
	peelVertices(level);//CtoCX
	if(control==1|| true){
		pruned=collectRemovedEV(level);
		if(!pruned){
			if(removeVE(level)) goto REC;
		}else goto REC;
	}
	if(MEInP>k || existSatV()|| CX_end<=LB){
		goto REC;
	}
	if(P_end>LB) {store(P_end);}
	if(MEInG <= k){store(CX_end); goto REC;}
#ifdef DPBOUND
	if(complexBound()<=LB) { goto REC;}
#endif
	tree_cnt++;
	int u;
	if (judgeDsetFull(PUB)){
		if(!isDSet()) goto REC;
		removeNonCN(level);
		maxClique4CN();
		goto REC;
	}
	u=chooseVertex();//pick a vertex u from C
	oldSz2=move_record.size();
	while (!Qe.empty()) Qe.pop();
	assert(Qv.empty());
	CtoP(u,level);
	if(false){
		collectRemovedEV(level);
		removeVE(level);
	}
	control=1;
	DSetBranch(level+1, PUB);
	if(false){
		recover(level, oldSz2);
	}
	PtoCX(u,level);//PtoCX(0,2)

	move_record.push_back(make_tuple(make_pair(-1, u), make_pair(0,2),-1));
	control=-1;
	DSetBranch(level+1, PUB);
	if(!isInCX(u)){
		printf("DSetBranch err, not in CX: subNum: %d, level: %d, treeIdx: %lld\n", subNum, level, treeIdx);
	}
	// CXtoC(u,level);//CXtoC
REC:
	if(useNew){
		recover(level, oldSz);
	}
}
void ExactSearcher::search(){
	int u=PC[0]; 
	szSum+=this->n;
	subNum++;
	double curDense=double(this->m/(n*(n-1)/2+0.0001));
	avgSubDense=((subNum-1)*avgSubDense+curDense)/subNum;
	collectRemovedEV(0);
	removeVE(0);
	if(isInX(0)) return;
	// CtoPwithPrune(u,0);
	CtoP(u,0);
	control=1;
	DSetBranch(1, 2*k+1);
	PtoC(u,0);
}

bool ExactSearcher::includeV(int level, int PUB){
	for (int i = P_end; i < C_end;i++ ){
		int u=PC[i];
		if(neiInG[u]==CX_end-2 && neiInP[u]==P_end-1 && MEInP<k) {
			include_level[u]= level;
			// restore_set[u]=1;
			CtoP(u, level);
			if(useNew){
				tuple<pair<int,int>,pair<int,int>,int> t(make_pair(-1,u), make_pair(1,0),-1);
				move_record.push_back(t);
			}
			PUB=min(PUB,P_end+2*(k-MEInP));
		}
		// else i++;
		if (judgeDsetFull(PUB)) return true;
	}
	
	return false;
}
void ExactSearcher::peelVertices(int level){
	for (int i = P_end; i < C_end; ){
		int u = PC[i];
		if(neiInP[u]+neiInC[u]==C_end-1 && peel_level[u]==n){
			// remove_level[u]=level;
			peel_level[u]=level;
			// restore_set[u]=1;
			if(useNew){
				tuple<pair<int,int>, pair<int,int>, int> t=make_tuple(make_pair(-1,u),make_pair(1,2),-1);
				// auto a=get<0>(t);
				move_record.push_back(t);
			}
			CtoCX(u, level);
		}
		else i++;
	}
}
bool ExactSearcher::collectRemovedEV(int level){
	bool check=true;
	if(check){
		//check the degree in P
		for (int i = 0; i < P_end; i++){
			int u=PC[i];
			if(P_end+(neiInG[u]-neiInP[u])+(k-MEInP)<=LB) return true;
		}
		
		//check the second order reduction
		if(enable){
			for (int i = 0; i < P_end; i++){
				for (int j = i+1; j < P_end; j++){
					if(edgePrune(PC[i],PC[j])) return true;
				}
			}
		}
	}
	

	for (int i = P_end; i < C_end;i++ ){
		int u = PC[i];
		if((P_end - neiInP[u] + MEInP > k || 1 + neiInG[u]+k-MEInP<=LB) && remove_level[u] == n){
			remove_level[u]=level;
			restore_set[u]=1;
			Qv.push(u);
			continue;
			// removed_v.push_back(u);
			// CtoX(u);
		}
		if(enable){
			// if(remove_level[u]<=level) continue;
			for (int j = 0; j < P_end; j++){
				int v=PC[j];
				if((edgePrune(v,u)|| MuEx[v][u]) && remove_level[u]==n){
					remove_level[u]=level;
					restore_set[u]=1;
					Qv.push(u);
					break;
				}
			}
		}
	}
	for (int i = C_end; i < CX_end;i++){
		int u = PC[i];
		if ((neiInP[u]<P_end || 1 + neiInG[u]+k-MEInP<=LB) && remove_level[u]==n){
			remove_level[u]=level;
			restore_set[u] = 2;
			Qv.push(u);
			continue;
			// removed_v.push_back(u);
			// CXtoX(u);
		}
		if(enable){
			// if(remove_level[u]<=level) continue;
			for (int j = 0; j < P_end; j++){
				int v=PC[j];
				if((edgePrune(v,u)|| MuEx[v][u]) && remove_level[u]==n){
					remove_level[u]=level;
					restore_set[u]=2;
					Qv.push(u);
					break;
				}
			}
		}
		// else i++;
	}

	//collect the mutual exclusive pairs and the delete edges
	int MuExNum=0;
	if(enable){
		for (int i = P_end; i < CX_end; i++){
			int u = PC[i];
			if(remove_level[u]<=level) continue;
			for (int j = i+1; j < CX_end; j++){
				int v=PC[j];
				if(remove_level[v]<=level) continue;
				if(edgePrune(u,v)) {
					MuExNum++;
					Qe.push(make_pair(u,v));
				}
			}
		}
	}
	double dense=(double)MuExNum/((C_end-P_end)*(C_end-P_end-1)/2+0.001);
	if(C_end-P_end>=5)denSum+=dense,denNum++;
	// denSum+=dense; denNum++;
	// printf("del density: %.2f%%\n",dense*100);
	MaxMuExNum=max(MaxMuExNum,MuExNum);
	MuExSum+=MuExNum;
	return false;
}
bool ExactSearcher::removeVE(int level){
	bool terminate=false;
	// remove the vertices in C
	// printf("enter removeVE\n");
	vector<int> u_neis;
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in removeVE begin\n", subNum, level);
		exit(0);
	}
#endif
	bool recurDelE=true;
	bool recurDel=true;
	if(useNew){
		//test
		long long ttsz=0;
		
		while (!Qv.empty() || !Qe.empty())
			{
				while (!Qv.empty()){
					int u=Qv.front();
					Qv.pop();
					if(restore_set[u]==1) {
						// CtoX(u);
						assert(isInC(u));
						swapID(PC_rid[u], --C_end);
						swapID(C_end, --CX_end);
						assert(isInX(u));
						assert(u==PC[CX_end]);
						int nonadj_inC = (C_end-P_end)-neiInC[u];
						int nonadj_inPC = C_end - (neiInP[u]+neiInC[u]);
						int nonadj_inG = CX_end - neiInG[u];
						int nonadj_inCCX = (CX_end - P_end) - (neiInG[u] - neiInP[u]);
						MEInC -= nonadj_inC;
						MEInPC -= nonadj_inPC;
						MEInG -= nonadj_inG;
						MEInCCX -= nonadj_inCCX;
						if(usemtx){
							vector<int> unei;
							unei.clear();
							for (int i = 0; i < CX_end; i++){
								int v=PC[i];
								if(isAdj(u,v)) unei.push_back(v);
							}
							for (int i = 0; i < unei.size(); i++){
								int v=unei[i];
								assert(isAdj(u,v));
								neiInC[v]--, neiInG[v]--;
								if(recurDel){
									if(isInP(v)){
										if(P_end+neiInG[v]-neiInP[v]+(k-MEInP)<=LB){
											terminate=true;
										}
									}
									else if(1 + neiInG[v]+k-MEInP<=LB && remove_level[v]==n){
										if(isInC(v)){
											Qv.push(v);
											remove_level[v]=level;
											restore_set[v]=1;
										}else if(isInCX(v)){
											Qv.push(v);
											remove_level[v]=level;
											restore_set[v]=2;
										}
									}
								}
							}
							if(terminate){
								// printf("terminate\n");
								swapID(PC_rid[u], CX_end++);
								swapID(CX_end-1, C_end++);
								remove_level[u]=n;
								MEInC += nonadj_inC;
								MEInPC += nonadj_inPC;
								MEInG += nonadj_inG;
								MEInCCX += nonadj_inCCX;
								for (int i = 0; i < unei.size(); i++){
									int v=unei[i];
									neiInC[v]++, neiInG[v]++;
								}
								// printf("u_0 neiInC: %d\n", neiInC[0]);
								goto CHECK;
								// return terminate;// return true;
							}
							if(enable){
								for (int id1 = 0; id1 < unei.size(); id1++)
								{
									int v=unei[id1];
									for (int id2 = id1+1; id2 < unei.size(); id2++)
									{
										int w=unei[id2];
										assert(isAdj(v,u) && isAdj(w, u));
		#ifdef _DBUG_
								int common_neighbors = 0;
					for(int k = P_end;k <= CX_end;k ++) if(matrix[PC[k]*n + v]&&matrix[PC[k]*n + w]) ++ common_neighbors;
					assert(CmNei[v][w] == common_neighbors);
					assert(CmNei[w][v] == common_neighbors);
		#endif
										CmNei[w][v]--;
										CmNei[v][w]--;
										if(recurDelE){
											if (edgePrune(v,w))
											{
												if(!isInP(v)&& !isInP(w)){
													// printf("enter recurdelE\n");
													Qe.push(make_pair(v,w));
												}else if(isInP(v) && !isInP(w)){
													if(remove_level[w]>level){
														remove_level[w]=level;
														Qv.push(w);
														if(isInC(w)) restore_set[w]=1;
														else restore_set[w]=2;
													}
												}else if(isInP(v) && isInP(w)){
													terminate=true;
												}
											}
											
										}
		#ifdef _DBUG_
					int common_neighbors2 = 0;
					for(int k = P_end;k < CX_end;k ++) if(matrix[PC[k]*n + v]&&matrix[PC[k]*n + w]) ++ common_neighbors2;
					if(CmNei[v][w]!=common_neighbors2){
						printf("subNum: %d, level: %d, u: %d, CmNei[%d][%d]=%d, common_neighbors2=%d, common_neighbors=%d\n",subNum,level,u, v,w,CmNei[v][w],common_neighbors2,common_neighbors);
					}
					assert(isAdj(v,u) && isAdj(w, u));
					assert(CmNei[v][w] == common_neighbors2);
					assert(CmNei[w][v] == common_neighbors2);
		#endif
									}
									
								}
							}
							unei.clear();
						}
						tuple<pair<int,int>, pair<int,int>, int> t(make_pair(-1, u), make_pair(1,3),-1);
						move_record.push_back(t);
					}
					else if(restore_set[u]==2){
						// CXtoX(u);
						assert(isInCX(u));
						swapID(PC_rid[u], --CX_end);
						assert(isInX(u));
						assert(u==PC[CX_end]);
						int nonadj_inG = CX_end - neiInG[u];
						int nonadj_inCCX = (CX_end - P_end) - (neiInG[u] - neiInP[u]);
						MEInG -= nonadj_inG;
						MEInCCX -= nonadj_inCCX;
						if (usemtx){
							vector<int> unei;
							for (int i = 0; i < CX_end; i++){
								int v=PC[i];
								if(isAdj(u,v)) unei.push_back(v);
							}
							for (int i = 0; i < unei.size(); i++){
								int v=unei[i];
								neiInG[v]--;
								if(recurDel){
									if(isInP(v)){
										if(P_end+neiInG[v]-neiInP[v]+(k-MEInP)<=LB){
											terminate=true;
										}
									}else if(1 + neiInG[v]+k-MEInP<=LB && remove_level[v]==n){
										if(isInC(v)){
											Qv.push(v);
											remove_level[v]=level;
											restore_set[v]=1;
										}else if(isInCX(v)){
											Qv.push(v);
											remove_level[v]=level;
											restore_set[v]=2;
										}
									}
								}
							}
							if(terminate){
								// printf("terminate\n");
								swapID(PC_rid[u], CX_end++);
								remove_level[u]=n;
								MEInG += nonadj_inG;
								MEInCCX += nonadj_inCCX;
								for (int i = 0; i < unei.size(); i++){
									int v=unei[i];
									neiInG[v]++;
								}
								goto CHECK;
								// return terminate;// return true;
							}
							if(enable){
								for (int id1 = 0; id1 < unei.size(); id1++)
								{
									int v=unei[id1];
									for (int id2 = id1+1; id2 < unei.size(); id2++)
									{
										int w=unei[id2];
		#ifdef _DBUG_
								int common_neighbors = 0;
					for(int k = P_end;k <= CX_end;k ++) if(matrix[PC[k]*n + v]&&matrix[PC[k]*n + w]) ++ common_neighbors;
					assert(CmNei[v][w] == common_neighbors);
					assert(CmNei[w][v] == common_neighbors);
		#endif
										CmNei[w][v]--;
										CmNei[v][w]--;
										if(recurDelE){
											if (edgePrune(v,w))
											{
												if(!isInP(v)&& !isInP(w)){
													// printf("enter recurdelE\n");
													Qe.push(make_pair(v,w));
												}else if(isInP(v) && !isInP(w)){
													if(remove_level[w]>level){
														remove_level[w]=level;
														Qv.push(w);
														if(isInC(w)) restore_set[w]=1;
														else restore_set[w]=2;
													}
												}else if(isInP(v) && isInP(w)){
													terminate=true;
												}
											}
											
										}
		#ifdef _DBUG_
					int common_neighbors2 = 0;
					for(int k = P_end;k < CX_end;k ++) if(matrix[PC[k]*n + v]&&matrix[PC[k]*n + w]) ++ common_neighbors2;
					if(CmNei[v][w]!=common_neighbors2){
						printf("subNum: %d, level: %d, CmNei[%d][%d]=%d, common_neighbors=%d\n",subNum,level, v,w,CmNei[v][w],common_neighbors2);
					}
					assert(CmNei[v][w] == common_neighbors2);
					assert(CmNei[w][v] == common_neighbors2);
		#endif
									}
									
								}
							}
							
							unei.clear();
						}
						tuple<pair<int,int>, pair<int,int>, int> t(make_pair(-1, u), make_pair(2,3),-1);
						move_record.push_back(t);
					}
					if(terminate){
						goto CHECK;
					}
				}
	
				if(Qe.empty()) break;

				int v1=Qe.front().first, v2=Qe.front().second;
				Qe.pop();
				if(!delEdge) continue;
				// printf("del edge\n");
				if(remove_level[v1]<=level|| remove_level[v2]<=level ||MuEx[v1][v2]) continue;
				// assert(isInC(v1));
				// assert(isInC(v2));
				MuEx[v1][v2]=MuEx[v2][v1]=true;
				int oriadj=isAdj(v1,v2)?1:0;
				if(isAdj(v1,v2)){
					isAdj(v1,v2)=isAdj(v2,v1)=false;
					neiInG[v1]--, neiInG[v2]--;
					MEInG++;
					if(isInC(v1) && isInC(v2)){
						neiInC[v1]--, neiInC[v2]--;
						MEInC++, MEInCCX++, MEInPC++;
					}else if(isInC(v1) && isInCX(v2)){
						MEInCCX++;neiInC[v2]--;
					}else if(isInCX(v1) && isInC(v2)){
						MEInCCX++;neiInC[v1]--;
					}else{
						MEInCCX++;
					}
					
					if(recurDel){
						// printf("enter recur del in edge del\n");
						if(1+neiInG[v2]+k-MEInP<=LB){
							if(isInC(v2)){
								restore_set[v2]=1;
								remove_level[v2]=level;
								Qv.push(v2);
							}
						}
						if(1+neiInG[v1]+k-MEInP<=LB){
							if(isInC(v1)){
								restore_set[v1]=1;
								remove_level[v1]=level;
								Qv.push(v1);
							}
						}
					}
					vector<int> unei;
					for (int id = 0; id < CX_end; id++){
						int v3=PC[id];
						if(isAdj(v1,v3)) unei.push_back(v3);
					}
					for (int i = 0; i < unei.size(); i++){
						int v3=unei[i];
						CmNei[v2][v3]--, CmNei[v3][v2]--;
						if(recurDelE){
							if(edgePrune(v2,v3)){
								if(!isInP(v2) && !isInP(v3)){
									Qe.push(make_pair(v2,v3));
								}else if(isInP(v2) && !isInP(v3)){
									if(remove_level[v3]>level){
										remove_level[v3]=level;
										Qv.push(v3);
										if(isInC(v3)) restore_set[v3]=1;
										else restore_set[v3]=2;
									}
								}else if(!isInP(v2) && isInP(v3)){
									if(remove_level[v2]>level){
										remove_level[v2]=level;
										Qv.push(v2);
										if(isInC(v2)) restore_set[v2]=1;
										else restore_set[v2]=2;
									}
								}
							}
						}
					}
					unei.clear();
					// vector<int> unei;
					for (int id = 0; id < CX_end; id++){
						int v3=PC[id];
						if(isAdj(v2,v3)) unei.push_back(v3);
					}
					for (int i = 0; i < unei.size(); i++){
						int v3=unei[i];
						CmNei[v1][v3]--, CmNei[v3][v1]--;
						if(recurDelE){
							if(edgePrune(v1,v3)){
								if(!isInP(v1) && !isInP(v3)){
									Qe.push(make_pair(v1,v3));
								}else if(isInP(v1) && !isInP(v3)){
									if(remove_level[v3]>level){
										remove_level[v3]=level;
										Qv.push(v3);
										if(isInC(v3)) restore_set[v3]=1;
										else restore_set[v3]=2;
									}
								}else if(!isInP(v1) && isInP(v3)){
									if(remove_level[v1]>level){
										remove_level[v1]=level;
										Qv.push(v1);
										if(isInC(v1)) restore_set[v1]=1;
										else restore_set[v1]=2;
									}
								}
							}
						}
					}
					unei.clear();
				}
				if(isInC(v1) && isInC(v2)){
					tuple<pair<int, int>,pair<int, int>,int> t(make_pair(v1,v2), make_pair(1,1), oriadj);
					move_record.push_back(t);
				}else if(isInC(v1) && isInCX(v2)){
					tuple<pair<int, int>,pair<int, int>,int> t(make_pair(v1,v2), make_pair(1,2), oriadj);
					move_record.push_back(t);
				}else if(isInCX(v1) && isInC(v2)){
					tuple<pair<int, int>,pair<int, int>,int> t(make_pair(v1,v2), make_pair(2,1), oriadj);
					move_record.push_back(t);
				}else{
					tuple<pair<int, int>,pair<int, int>,int> t(make_pair(v1,v2), make_pair(2,2), oriadj);
					move_record.push_back(t);
				}
			}
			
	}

CHECK:
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in removeVE end\n", subNum, level);
		exit(0);
	}
#endif
	u_neis.clear();
	return terminate;
}
void ExactSearcher::recover(int level, int oldSz){
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in recover begin\n", subNum, level);
		exit(0);
	}
#endif
// printf("recover new\n");
	while(!Qv.empty()){
		int u=Qv.front();
		remove_level[u]=n;
		Qv.pop();
	}
	while (!Qe.empty()) Qe.pop();
	// while (!Qex.empty()) Qex.pop();
	int popSz=move_record.size()- oldSz;
	assert(popSz>=0);
	while (popSz>0){
		popSz--;
		tuple<pair<int,int>, pair<int,int>, int> t=move_record.back();
		auto e=get<0>(t);
		
		if(e.first==-1){
			int u=e.second;
			auto track=get<1>(t);
			if(track.first==1 && track.second==2){
				peel_level[u]=n;
				CXtoC(u,level);
			}else if(track.first==1 && track.second==0){
				include_level[u]=n;
				PtoC(u,level);
			}else if(track.first==1 && track.second==3){
				remove_level[u]=n;
				XtoC(u,level);
			}else if(track.first==0 && track.second==2){
				CXtoC(u,level);
			}else{
				remove_level[u]=n;
				XtoCX(u,level);
			}
		}else{
			if(!delEdge){
				//if not delete the edges
				int u=e.first,v = e.second;
				MuEx[u][v]=MuEx[v][u]=false;
			}else{
				int u=e.first,v=e.second;
				auto adj=get<2>(t);
				MuEx[u][v]=MuEx[v][u]=false;
				if(adj==1){
					auto loc=get<1>(t);
					isAdj(u,v)=isAdj(v,u)=true;
					neiInG[u]++, neiInG[v]++;
					MEInG--;

					if(loc.first==1 && loc.second==1){
						//both in C
						neiInC[u]++, neiInC[v]++;
						MEInC--, MEInCCX--, MEInPC--;
					}else if(loc.first==1 && loc.second==2){
						MEInCCX--;neiInC[v]++;
					}else if(loc.first==2 && loc.second==1){
						MEInCCX--;neiInC[u]++;
					}else{
						MEInCCX--;
					}
					
					for (int i = 0; i < CX_end; i++){
						int w=PC[i];
						if(isAdj(v,w)) CmNei[u][w]++, CmNei[w][u]++;
					}
					for (int i = 0; i < CX_end; i++){
						int w=PC[i];
						if(isAdj(u,w)) CmNei[v][w]++, CmNei[w][v]++;
					}
				}
			}
		}
		move_record.pop_back();
	}
#ifdef _CHECK_COIN_
	if(checkConinfidence(level)==-1){
		printf("subNum: %d, level: %d, error in recover end\n", subNum, level);
		exit(0);
	}
#endif

	#ifdef _DBUG_
	vector<int> P_rcv;
	vector<int> C_rcv;
	vector<int> CX_rcv;

	int MEInPRcv=MEInP;
	int MEInCRcv=MEInC;
	int MEInGRcv=MEInG;
	int MEInPCRcv=MEInPC;
	int MEInCCXRcv=MEInCCX;

	vecPush(P_rcv,PC,0,P_end);
	vecPush(C_rcv,PC,P_end,C_end);
	vecPush(CX_rcv,PC,C_end, CX_end);
	if(!vecCmp(P_rcd[level], P_rcv)){
		printf("subNum: %d, level: %d, treeIdx: %lld, P incoincidence\n",subNum, level ,treeIdx);
		exit(0);
	}
	assert(vecCmp(P_rcd[level], P_rcv));
	if(!vecCmp(C_rcd[level], C_rcv)){
		printf("subNum: %d, level: %d, treeIdx: %lld, C incoincidence\n",subNum, level,treeIdx);
		exit(0);
	}
	assert(vecCmp(C_rcd[level], C_rcv));
	if(!vecCmp(CX_rcd[level], CX_rcv)){
		printf("subNum: %d, level: %d, treeIdx: %lld, CX inconincidence\n",subNum, level, treeIdx);
		exit(0);
	}
	assert(vecCmp(CX_rcd[level], CX_rcv));

	if(MEInPRcd[level]!=MEInPRcv){
		printf("recover err, MEInP, subNum: %d, level: %d, treeIdx: %d\n", subNum, level, treeIdx);
		exit(0);
	}
	if(MEInCRcd[level]!=MEInCRcv){
		printf("recover err, MEInC, subNum: %d, level: %d, treeIdx: %d\n", subNum, level, treeIdx);
		exit(0);
	}
	if(MEInPCRcd[level]!=MEInPCRcv){
		printf("recover err, MEInPC, subNum: %d, level: %d, treeIdx: %d\n", subNum, level, treeIdx);
		exit(0);
	}
	if(MEInCCXRcd[level]!=MEInCCXRcv){
		printf("recover err, MEInCCX, subNum: %d, level: %d, treeIdx: %d\n", subNum, level, treeIdx);
		exit(0);
	}
	if(MEInGRcd[level]!=MEInGRcv){
		printf("recover err, MEInG, subNum: %d, level: %d, treeIdx: %d\n", subNum, level, treeIdx);
		exit(0);
	}

	
	P_rcv.clear();
	C_rcv.clear();
 	CX_rcv.clear();

#endif
	// printf("recover new end\n");
}
void ExactSearcher::removeNonCN(int level){
	// printf("enter remove\n");
	for (int i = P_end; i < C_end; ){
		int u = PC[i];
		if (neiInP[u]<P_end && remove_level[u]==n){
			remove_level[u]=level;
			CtoX(u, level);
			if (useNew){
				tuple<pair<int,int>, pair<int,int>, int > t(make_pair(-1, u), make_pair(1,3),-1);
				move_record.push_back(t);
			}
			
		}
		else i++;
	}
	
	for (int i = C_end; i < CX_end;){
		int u = PC[i];
		if (neiInP[u]<P_end && remove_level[u]==n){
			remove_level[u]=level;
			CXtoX(u, level);
			if(useNew){
				tuple<pair<int,int>, pair<int,int>, int > t(make_pair(-1, u), make_pair(2,3),-1);
				move_record.push_back(t);
			}
		}
		else i++;
	}
}

bool ExactSearcher::judgeDsetFull(int P_UB){
	// printf("enter judgedestfull\n");
	bool flag=false;
	if (MEInP==k || (MEInPC==MEInP) ||P_end==C_end ||P_end>=P_UB||P_end == 2*k+1 ){//P_end==C_end can be deleted
		flag = true;
	}
	// printf("exit judgedsetfull\n");
	return flag;
}
int ExactSearcher::bound(int neis){
	return min(P_end+neis+(k-MEInP),UB);
}
int ExactSearcher::candiBucket(){
	int CMax=0, kCur=k-MEInP;
	vector<vector<int>> nonNeiBucket(kCur+1);
	for (int i = P_end; i < CX_end; i++){
		int u=PC[i];
		if(P_end-neiInG[u]>kCur) continue;
		nonNeiBucket[P_end-neiInP[u]].push_back(u);
	}
	for (int i = 0; i <= kCur; i++){
		if(kCur<i) goto RET;
		for (auto u:nonNeiBucket[i]){
			if(kCur<i) goto RET;
			kCur-=i;
			CMax++;
		}
	}
RET:
	nonNeiBucket.clear();
	return P_end+CMax;
}
int ExactSearcher::PackingBound(){
	int CMax=0, kCur=k-MEInP;
	vector<int> nonNeiBucket(kCur+1);
	for (int i = P_end; i < CX_end; i++){
		int u=PC[i];
		if(P_end-neiInG[u]>kCur) continue;
		nonNeiBucket[P_end-neiInP[u]]++;
	}
	for (int i = 0; i <= kCur; i++){
		if(kCur-i*nonNeiBucket[i]>=0) CMax+=nonNeiBucket[i], kCur-=i*nonNeiBucket[i];
		else{CMax+=kCur/i;break;}
	}
RET:
	nonNeiBucket.clear();
	return P_end+CMax;
}
int ExactSearcher::ColorBound(){
		int C_max = 0, k_cur=k-MEInP; int maxCol=0;
		int i; int cnt=0;
		colorSz[0]=0;int col=0;
		//color CX first
		for (int i = C_end; i < CX_end; i++){
			int u=PC[i];
			if(P_end>neiInP[u]) continue;
			col=0;
			while (colorUseMtx[u*n+col]==treeIdx) col++;
			if(col>maxCol) maxCol=col, colorSz[maxCol]=0;
			colVecC[n*col+colorSz[col]]=u;
			colorSz[col]++;
			//maintain the color use
			for (int j = P_end; j < CX_end; j++){
				if(isAdj(u,PC[j])) colorUseMtx[n*PC[j]+col]=treeIdx;
			}
			colorLabel[u]=col;
		}

		//color C second
		for (int i = P_end; i < C_end; i++){
			col=0;
			int u=PC[i];
			while (colorUseMtx[u*n+col]==treeIdx) col++;
			if(col>maxCol) maxCol=col, colorSz[maxCol]=0;
			colVecC[n*col+colorSz[col]]=u;
			colorSz[col]++;
			for (int j = P_end; j < C_end; j++){
				if(isAdj(u,PC[j])) colorUseMtx[n*PC[j]+col]=treeIdx;
			}
			colorLabel[u]=col;
		}
		
	// colorNum=maxCol+1;

		while( k_cur>= cnt) {
			for(int j=0;j<=maxCol;j++){
				if(colorSz[j]>0 && k_cur>=cnt){
					k_cur-=cnt;
					--colorSz[j];
					++C_max;
				}
			}
			++cnt;
		}
		return P_end+C_max;

}
int ExactSearcher::CLUB(){
	int MaxNum=0, colorNum=0, maxNonDeg=0, maxWeight=0;//MaxNum- max num of chosen from C; 
	int kCur=k-MEInP, maxCol=0, col=0;
	colorSz[0]=0;
	for (int i = P_end; i < CX_end; i++){
		int u=PC[i];
		if(P_end-neiInP[u]>kCur) continue;
		// assert(P_end-neiInP[u]>0);
		nonNeiInPB[P_end-neiInP[u]].push_back(u);
		maxNonDeg=max(maxNonDeg, P_end-neiInP[u]);
	}
	for (int i = 0; i <= maxNonDeg; i++){
		if(kCur-i*nonNeiInPB[i].size()>=0) MaxNum+=nonNeiInPB[i].size(), kCur-=i*nonNeiInPB[i].size();
		else{MaxNum+=kCur/i;break;}
	}
	if(P_end+MaxNum<=LB) goto CLEAR;

	MaxNum=0, kCur=k-MEInP;
	for (int i = 0; i <= maxNonDeg; i++){
		int remain=nonNeiInPB[i].size(), initW=i;
		for (int j = 0; j < nonNeiInPB[i].size(); j++){
			int u=nonNeiInPB[i][j];bool flag=false;
			col=0;
			while(col<=maxCol){
				for (int id = 0; id < colorSz[col]; id++){
					int v=colVecC[n*col+id];
					if(isAdj(u,v)){col++,flag=true;break;}
				}
				if(!flag) break;
				flag=false;
			}
			if(col>maxCol){
				maxCol=col;
				colorSz[maxCol]=0;
			}
			colVecC[n*col+colorSz[col]]=u;
			colorSz[col]++;
		}
		colorNum=maxCol+1;
		//put the weight into the weight-bucket
		while(remain>0){
			weightBucket[initW]+=min(colorNum, remain);
			remain-=colorNum;
			maxWeight=max(maxWeight,initW);
			initW++;
		}
		maxCol=0, colorSz[0]=0,colorNum=0;//reset
	}
	for (int w = 0; w <= maxWeight; w++){
		if(kCur-w*weightBucket[w]>=0) MaxNum+=weightBucket[w], kCur-=w*weightBucket[w];
		else {
			MaxNum+=kCur/w;
			break;
		}
	}
	
CLEAR:
	for (int d = 0; d<= maxNonDeg; d++) nonNeiInPB[d].clear();// clear the nonNeiInPB
	memset(weightBucket,0, (maxWeight+1)*sizeof(int));
	return P_end+MaxNum;
}
int ExactSearcher::complexBound(){
	int MaxNum=0, colorNum=0, maxNonDeg=0, maxWeight=0;//MaxNum- max num of chosen from C; 
	totExNum=0;
	int kCur=k-MEInP, maxCol=-1, col=0;
	bool useGroupColor=false;
	//use bucket to collect the vertices in C
	for (int i = P_end; i < C_end; i++){
		int u=PC[i];
		if(P_end-neiInP[u]>kCur) continue;
		// assert(P_end-neiInP[u]>0);
		nonNeiInPB[P_end-neiInP[u]].push_back(u);
		maxNonDeg=max(maxNonDeg, P_end-neiInP[u]);
	}
	
	//color CX first
	for (int i = C_end; i < CX_end; i++){
		int u=PC[i];
		if(P_end>neiInP[u]) continue;
		col=0;
		while (colorUseMtx[u*n+col]==treeIdx) col++;
		if(col>maxCol) maxCol=col, colorSz[maxCol]=0, colorSzCX[maxCol]=0, colExNum[maxCol]=0;
		colVecCX[n*col + colorSzCX[col]]=u;
		colorSzCX[col]++;
		totExNum+=colorSzCX[col]-1;// sum the exclusive pairs
		//maintain the color use
		for (int j = P_end; j < CX_end; j++){
			if(isAdj(u,PC[j])) colorUseMtx[n*PC[j]+col]=treeIdx;
		}
		colorLabel[u]=col;
	}
	
	//then color C
	for (int d = 0; d <= maxNonDeg; d++){
		for(auto u:nonNeiInPB[d]){
			col=0;
			while (colorUseMtx[u*n+col]==treeIdx) col++;
			if(col>maxCol) maxCol=col, colorSz[maxCol]=0, colorSzCX[maxCol]=0, colExNum[maxCol]=0;
			colVecC[n*col+colorSz[col]]=u;
			colorSz[col]++;
			totExNum+=colorSzCX[col];
			for (int j = 0; j < colorSz[col]-1; j++){
				if(MuEx[u][colVecC[n*col+j]]) totExNum++, colExNum[col]++;
			}
			for (int j = P_end; j < C_end; j++){
				if(isAdj(u,PC[j])) colorUseMtx[n*PC[j]+col]=treeIdx;
			}
			colorLabel[u]=col;
		}
	}
	colorNum=maxCol+1;
	for (int c = 0; c <= maxCol; c++){
		//first assign weight for the CX
		for (int i = 0; i < colorSzCX[c]; i++){
			int u=colVecCX[n*c+i];
			weight[u]=P_end-neiInP[u]+i;
			maxWeight=max(weight[u],maxWeight);
			weightBucket[weight[u]]++;
		}
		//then assign weight for the C
		for (int i = 0; i < colorSz[c]; i++){
			int u=colVecC[n*c+i];
			weight[u]=P_end-neiInP[u]+colorSzCX[c]+i;
			maxWeight=max(maxWeight, weight[u]);
			weightBucket[weight[u]]++;
		}
	}
	for (int w = 0; w <= maxWeight; w++){
		if(kCur-w*weightBucket[w]>=0) MaxNum+=weightBucket[w], kCur-=w*weightBucket[w];
		else {
			MaxNum+=kCur/w;
			break;
		}
	}
	if(P_end+MaxNum<=LB) goto CLEAR;

	if(totExNum>0){
		kCur=k-MEInP;
		//init t(c,x)
		for(int c = 0; c<=maxCol; c++){
			int sz=0,ws=0;//ws-weight sum
			if(colExNum[c]>0) useGroupColor=true;
			if(useGroupColor){
				int* colvec=colVecC+n*c;
				colPacker->set(MuEx, colvec, colorSz[c]);
				colPacker->coloring(colvec,neiInP, P_end);
				for(int k1=0; k1<=kCur; k1++){
					t[c*n+k1]=colorSzCX[c]>0?1:0;
					//compute the sz
					sz=colPacker->colPack(colvec, colorSz[c], neiInP, P_end, k1);
					t[c*n+k1]=max(t[c*n+k1],sz);
					ws=0, sz=0;
				}
				colPacker->reset();
				sz=0;
			}else{
				for(int k1=0; k1<=kCur; k1++){
					t[c*n+k1]=colorSzCX[c]>0?1:0;
					//compute the sz
					for (int i = 0; i < colorSz[c]; i++){
						int u=colVecC[c*n+i];
						ws+=i+P_end-neiInP[u];
						if(ws<=k1) sz++;
					}
					t[c*n+k1]=max(t[c*n+k1],sz);
					ws=0, sz=0;
				}
			}
			useGroupColor=false;
		}

		// init DP(0,*)
		for (int i = 0; i <= kCur; i++){
			dp[0*n+i]=t[0*n+i];
			// DP(0,i)=T(0,i);
		}
		//do the dynamic programming
		for (int c = 1; c <= maxCol; c++){
			for (int k1 = 0; k1 <= kCur; k1++){	
				dp[c*n+k1]=0;
				for (int k2 = 0; k2 <= k1; k2++){
					dp[c*n+k1]=max(dp[c*n+k1],dp[(c-1)*n+k2]+t[c*n+(k1-k2)]);
					if(dp[c*n+k1]>LB-P_end){MaxNum=n;goto CLEAR;}
				}
			}
		}
		MaxNum=min(MaxNum, dp[maxCol*n+kCur]);
	}
CLEAR:
	for (int d = 0; d<= maxNonDeg; d++) nonNeiInPB[d].clear();// clear the nonNeiInPB
	memset(weightBucket,0, (maxWeight+1)*sizeof(int));
	return P_end+MaxNum;
}
bool ExactSearcher::bound(){//a simple upper bound
	int neis=0;
	for (int i = P_end; i < CX_end; i++){
		int v=PC[i];
		if(neiInP[v]==P_end) neis++;
	}
	return P_end+neis+(k-MEInP)<=LB;
}

void ExactSearcher::PrintColoring(){
	//print the color num
	printf("C_num=%d, CX_num=%d, colorNum=%d\n",C_end-P_end, CX_end-C_end,colNum);
	//print the coloring of C
	printf("color in C\n");
	int C_verify=0, CX_verify=0;
	for (int c = 0; c < colNum; c++){
		printf("color %d, colornum=%d: ", c, colorSz[c]);
		C_verify+=colorSz[c];
		for (int i = 0; i < colorSz[c]; i++){
			int u = colVecC[c*n+i];
			printf("%d %d %d, ", u, P_end-neiInP[u], colB[n*c+i]);
		}
		printf("\n");
	}
	//print the coloring of CX
	printf("color in CX\n");
	for (int c = 0; c < colNum; c++){
		CX_verify+=colorSzCX[c];
		printf("color %d, colornum=%d, colCX=%d: ", c, colorSzCX[c], colCX[c]);
		for (int i = 0; i < colorSzCX[c]; i++){
			int u = colVecCX[c*n+i];
			printf("%d %d %d, ", u, P_end-neiInP[u], colB[n*c+i]);
		}
		printf("\n");
	}
	// assert(C_verify==C_end-P_end);
	// assert(CX_verify==CX_end-C_end);
}
void ExactSearcher::store(int newLB){
	KDC.resize(LB=newLB);
	maxME=max(maxME, MEInP);
	for (int i = 0; i < LB; i++) KDC[i]=PC[i];
}
void ExactSearcher::maxClique4CN(){
	subSearcher->loadCN(this);
	subSearcher->search();
	if(subSearcher->MC.size()+P_end>LB){
		KDC.resize(LB=subSearcher->MC.size()+P_end);
		for (int i = 0; i < P_end; i++) KDC[i]=PC[i];
		for (int i = 0; i < subSearcher->MC.size(); i++) KDC[i+P_end]=PC[subSearcher->MC[i]+P_end];
	}
}
bool ExactSearcher::vertexCmp(int u, int v){
	if(neiInP[u]<neiInP[v]) return true;
	if(neiInP[u]==neiInP[v] && neiInC[u]<neiInC[v]) return true;
	return false;
}
int ExactSearcher::chooseVertex(){
	int u = PC[P_end];
	for(int id=P_end;id<C_end;id++){
		int v=PC[id];
		if(vertexCmp(v,u))u=v;
	}
	return u;
}
void ExactSearcher::degreeSumPrune(int level){
	int kCur=k-MEInP;int degSum = 0;
	int id = P_end;
	vector<vector<int>> nonNeiBucket(kCur+1);
	for (int i = P_end; i < CX_end; i++)
	{
		int u=PC[i];
		nonNeiBucket[P_end-neiInP[u]].push_back(u);
	}

	for (int i = 0; i < nonNeiBucket.size(); i++){
		for(auto u: nonNeiBucket[i]){
			if(id < LB) degSum+=i, id++;
			else if(i > kCur - degSum){
				remove_level[u]=level;
				restore_set[u]=1;
				CtoX(u,level);
			}
		}
	}
	nonNeiBucket.clear();
	return;
}
void ExactSearcher::colorBasedPrune(int level){
	int weiSum = MEInP;
	int id = P_end;
	int maxCol=0;int maxW = 0;
	int kCur=k-MEInP;
	vector<vector<int>> nonNeiBucket(kCur+1);
	for (int i = P_end; i < CX_end; i++){
		int u = PC[i]; int c = colorLabel[u];
		maxCol=max(maxCol, c);
		if(P_end-neiInP[u]>kCur) continue;
		nonNeiBucket[P_end-neiInP[u]].push_back(u);
	}
	vector<vector<int>> weiBucket(kCur+n);
	vector<vector<int>> colBucket(maxCol+1); 
	for (int i = 0; i < nonNeiBucket.size(); i++){
		for(auto u: nonNeiBucket[i]){
			int c = colorLabel[u];
			int wei = P_end-neiInP[u]+colBucket[c].size();
			maxW=max(maxW, wei);
			colBucket[c].push_back(u), weiBucket[wei].push_back(u);
		}
	}
	for (int w = 0; w <= maxW ; w++){
		for(auto u:weiBucket[w]){
			if(id < LB) weiSum+=w, id++;
			else if(P_end-neiInP[u] > k-weiSum && isInC(u)){
				remove_level[u]=level;
				CtoX(u,level);
				if(useNew){
					tuple<pair<int,int>, pair<int,int>,int> t(make_pair(-1,u),make_pair(1,3),-1);
					move_record.push_back(t);
				}
			}
		}
	}
	nonNeiBucket.clear();
	colBucket.clear();
	weiBucket.clear();
	return;
}

void ExactSearcher::printPos(int u){
	printf("P_end: %d, C_end: %d, CX_end: %d, X_end: %d\n", P_end, C_end, CX_end, n);
	if (isInP(u)) printf("u: %d is in P\n");
	else if(isInC(u)) printf("u: %d is in C\n");
	else if(isInCX(u)) printf("u: %d is in CX\n");
	else printf("u: %d is in X\n");
}
bool ExactSearcher::existSatV(){
	for (int i = 0; i < P_end; i++){//scan the vertex in P
		int u = PC[i];
		if (neiInP[u]==C_end-1 && u!=0) return true;
	}
	return false;
}
bool ExactSearcher::isDSet(){
	for (int i = 0; i < P_end; i++){//scan the vertex in P
		int u = PC[i];
		if (neiInP[u]==P_end-1 && u!=0) return false;
	}
	return true;
}
void ExactSearcher::checkColUseMtx(){
	for (int i = P_end; i < C_end; i++){
		int u=PC[i];
		for (int c = 0; c < n; c++){
			if(colorUseMtx[u*n+c]>=treeIdx){
				printf("error of colUseMtx: %d, %d\n", subNum, treeIdx);
				assert(colorUseMtx[u*n+c]<treeIdx);
				exit(0);
			}
		}
	}
}

void ExactSearcher::swapID(int i, int j) {
	swap(PC[i], PC[j]);
	PC_rid[PC[i]] = i;
	PC_rid[PC[j]] = j;
}

SubSearcher::SubSearcher() {
	n = m = 0;

	pstart = nullptr;
	edges = nullptr;
	pend = nullptr;

	pstart_o = nullptr;
	edges_o = nullptr;

	tmp_vs = nullptr;
	tmp_color = nullptr;

	matrix = nullptr;

	head = nullptr;
	next = nullptr;

	maxCore = 0;

	vis = nullptr;

	degree = nullptr;
	rid = nullptr;
	id = nullptr;

	current_clique = nullptr;

	mapping = nullptr;
	mapping_n = matrix_len = 0;

	branches = 0;
}

SubSearcher::~SubSearcher() {
	release();
}

void SubSearcher::load(ExactSearcher* exactSearcher){
	const int* PC=exactSearcher->PC;
	const int P_end = exactSearcher->P_end;
	const int C_end = exactSearcher->C_end;
	const int EInC=(C_end-P_end)*(C_end-P_end-1)-2*(exactSearcher->MEInC);

	n=C_end-P_end, m= EInC;
	pstart = new ept[n+1]; edges = new ui[m]; pstart[0] = 0;
	
	int idx=0;
	for(int i=0;i<n;++i){
		pstart[i]=idx;	
		for(int j=0;j<n;++j){
			if(exactSearcher->isAdj(PC[P_end+i],PC[P_end+j])) edges[idx++]=j;
		}
	}
	pstart[n]=idx;

	MC.resize(exactSearcher->LB-P_end);
	// printf("\tn = %d; m = %ld (undirected)\n", n, m/2);
}

void SubSearcher::loadCN(ExactSearcher* exactSearcher){
	int pstartSz=sizeof(pstart)/sizeof(ept), edgeSz=sizeof(edges)/sizeof(ui);
	const int* PC = exactSearcher->PC;
	const int* neiInP = exactSearcher->neiInP;
	const int P_end = exactSearcher->P_end;
	const int CX_end = exactSearcher->CX_end;
	const int EInCCX = (CX_end - P_end) * (CX_end - P_end - 1) - 2 * (exactSearcher->MEInCCX);
	int idx = 0;
	this->n = CX_end - P_end; this->m = EInCCX;
	pstart = new ept[n+1]; 
	if (m<0){
		printf("P_end: %d, CX_end: %d, MEInCCX: %d\n",P_end, CX_end, exactSearcher->MEInCCX);
	}
	edges = new ui[m]; 
	// printf("new edges\n");
	pstart[0] = 0;
	for (int i = 0; i < n; i++){
		int u = PC[P_end + i];
		pstart[i] = idx;
		if (neiInP[u] < P_end) continue;
		for (int j = 0; j < n; j++){
			int v = PC[P_end + j];
			if(neiInP[v] == P_end && exactSearcher->isAdj(u,v)) edges[idx++]=j;
		}
		
	}
	pstart[n]=idx;
	this->m=idx;
	MC.resize(exactSearcher->LB-P_end);
	
}

void SubSearcher::read(const char *_file) {
	FILE *f = fopen(_file, "rb"); 
	ui tt; fread(&tt, sizeof(int), 1, f); fread(&n, sizeof(int), 1, f); fread(&m, sizeof(int), 1, f);
	ui *degree = new ui[n]; fread(degree, sizeof(int), n, f); pstart = new ept[n+1]; edges = new ui[m]; pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {
		fread(edges+pstart[i], sizeof(int), degree[i], f);
		ui *buff = edges+pstart[i]; sort(buff, buff+degree[i]); ui idx = 0;
		for(ui j = 0;j < degree[i];j ++) {
			if(buff[j] == i||(j > 0&&buff[j] == buff[j-1])) continue;
			buff[idx ++] = buff[j];
		}
		pstart[i+1] = pstart[i] + idx;
	}
	delete[] degree; fclose(f);
	// printf("\tn = %d; m = %ld (undirected)\n", n, m/2);
}

void SubSearcher::write() {
	FILE *fout = fopen("clique.txt", "w"); fprintf(fout, "%lu\n", MC.size());
	sort(MC.begin(), MC.end()); for(auto ele : MC) fprintf(fout, "%u ", ele);
	fclose(fout);
}

void SubSearcher::search() {
	// Timer t;
	ui *seq = new ui[n];
	ui *core = new ui[n];
	ui *color = new ui[n];
	ui UB = degen(seq, core, color);
	ui oldSz = MC.size();
	if(MC.size() < UB) {
		ui *out_mapping  = new ui[n];
		shrink_graph(seq, core, color, out_mapping, pstart_o, edges_o);
		current_clique = new ui[UB];
		mapping = new ui[maxCore];
		head = new ui[maxCore];
		next = new ui[maxCore];
		id = new ui[maxCore];
#ifdef _BITSET_
		matrix = new unsigned char[maxCore*((maxCore+7)/8)];
#else
		matrix = new unsigned char[max_core*max_core];
#endif
		degree = new ui[n];
		rid = new ui[n];
		vis = new char[n];
		memset(vis, 0, sizeof(char)*n);

		ui *subUB = new ui[n];
		UB = ego_degen(seq, core, color, subUB, UB);
		// printf("\tMC-EGO Time: %s\n", Utility::integer_to_string(t.elapsed()).c_str());

		if(UB > MC.size()) search_oriented(seq, core, color, subUB);

		if(MC.size() > oldSz) for(ui i = 0;i < MC.size();i ++) MC[i] = out_mapping[MC[i]];

		delete[] out_mapping;
		delete[] subUB;
	}

	delete[] color;
	delete[] core;
	delete[] seq;
	release();

	// printf("\tMaximum Clique Size: %lu, Max Depth: %u, Total Time: %s\n", MC.size(), max_depth, Utility::integer_to_string(t.elapsed()).c_str());
}

ui SubSearcher::coloring_csr(const ui *vs, const ui vs_size, const ui original_size, ui *color, char *vis, const ui start_idx, const ui start_color) {
	for(ui i = 0;i < original_size;i ++) color[vs[i]] = n;
	for(ui i = vs_size - start_color;i < vs_size;i ++) color[vs[i]] = vs_size - i - 1;

	ui max_color = 0;
	for(ui i = vs_size - start_color;i > start_idx;i --) {
		ui u = vs[i-1];
		for(ui j = pstart[u];j < pend[u];j ++) {
			ui c = color[edges[j]];
			if(c != n) vis[c] = true;
		}
		for(ui j = 0;;j ++) if(!vis[j]) {
			color[u] = j;
			if(j > max_color) max_color = j;
			break;
		}
		for(ui j = pstart[u];j < pend[u];j ++) {
			ui c = color[edges[j]];
			if(c != n) vis[c] = false;
		}
	}

	return max_color + 1;
}

ui SubSearcher::degen(ui *seq, ui *core, ui *color) {
	// heuristic_max_clique_max_degree(10);
	// Timer t;
	ui threshold = MC.size(); ui *degree = new ui[n];
	char *vis = new char[n]; memset(vis, 0, sizeof(char)*n);
	pend = new ept[n]; for(ui i = 0;i < n;i ++) pend[i] = pstart[i+1], degree[i] = pend[i] - pstart[i];

	ui queue_n = 0, new_size = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] < threshold) seq[queue_n ++] = i;
	for(ui i = 0;i < queue_n;i ++) {
		ui u = seq[i]; degree[u] = 0;
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(degree[edges[j]] > 0) {
			if((degree[edges[j]] --) == threshold) seq[queue_n ++] = edges[j];
		}
	}
	for(ui i = 0;i < n;i ++) {
		if(degree[i] >= threshold) seq[queue_n + (new_size ++)] = i;
		else vis[i] = 1, core[i] = 0;
	}

	ui UB = n;
	if(new_size == 0) UB = MC.size();
	else {
		ListLinearHeap *heap = new ListLinearHeap(n, new_size-1);
		heap->init(new_size, new_size-1, seq+queue_n, degree);
		maxCore = 0;
		vector<ui> res;
		for(ui i = 0;i < new_size;i ++) {
			ui u, key; heap->pop_min(u, key);
			if(key > maxCore) maxCore = key; core[u] = maxCore; seq[queue_n + i] = u;
			if(key + i + 1 == new_size) {
				ui x_size = i+1;
				heap->get_ids(seq+queue_n, x_size);
				assert(x_size == new_size);
				for(ui j = i;j < new_size;j ++) {
					core[seq[queue_n+j]] = maxCore;
					res.pb(seq[queue_n+j]);
				}
				break;
			}
			vis[u] = true;
			for(ui j = pstart[u];j < pstart[u+1];j ++) if(!vis[edges[j]]) heap->decrement(edges[j], 1);
		}
		delete heap;
		// printf("*** Degeneracy clique size: %lu, max_core: %u, Time: %s (microseconds)\n", res.size(), maxCore, Utility::integer_to_string(t.elapsed()).c_str());

		if(res.size() > MC.size()) MC = res;
		if(MC.size() == maxCore+1) UB = maxCore + 1;
		else{
			memset(vis, 0, sizeof(char)*n);
			ui start_idx = 0;
			while(start_idx < n && core[seq[start_idx]] < MC.size()) ++ start_idx;
			UB = coloring_csr(seq, n, n, color, vis, start_idx, 0);
			if(MC.size() < UB) {
				for(ui i = 0;i < res.size();i ++) vis[res[i]] = true;
				for(ui i = n-res.size();i > 0;i --) {
					ui u = seq[i-1], cnt = 0;
					if(core[u] < res.size()) break;
					for(ui j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]]) ++ cnt;
					if(cnt == res.size()) res.pb(u), vis[u] = true;
				}
				// printf("*** Degen_greedy_extend clique size: %lu, UB: %u, Time: %s (microseconds)\n", res.size(), UB, Utility::integer_to_string(t.elapsed()).c_str());
			}
			if(res.size() > MC.size()) MC = res;
		}
	}
	delete[] degree;
	delete[] vis;
	return UB;
}

ui SubSearcher::degeneracy_maximal_clique_adj_list(const ui current_clique_size, ui *vs, const ui vs_size, char *vis, ui *degree, ListLinearHeap *heap) {
	assert(vs_size > 0);

	for(ui j = 0;j < vs_size;j ++) vis[vs[j]] = 1;

	char sparse = 0;
	ui total_edges = 0;
	for(ui i = 0;i < vs_size;i ++) total_edges += degree[vs[i]];
	if(total_edges*10 < vs_size*(vs_size-1)) sparse = 1;

	if(sparse) heap->init(vs_size, vs_size-1, vs, degree);

	ui start_color = 0;
	for(ui j = 0;j < vs_size;j ++) {
		ui v, key;
		if(sparse) {
			heap->pop_min(v, key);
			vs[j] = v;
		}
		else {
			ui min_idx = j;
			for(ui k = j+1;k < vs_size;k ++) if(degree[vs[k]] < degree[vs[min_idx]]) min_idx = k;
			if(min_idx != j) swap(vs[min_idx], vs[j]);
			v = vs[j]; key = degree[v];
		}
		if(key + j + 1 == vs_size) {
			start_color = key + 1;
			if(sparse) {
				ui new_size = j+1;
				heap->get_ids(vs, new_size);
				assert(new_size == vs_size);
			}
			if(key + 1 + current_clique_size > MC.size()) {
				//printf("Find clique of size %u after search %u egos\n", degree[v] + 2, n-i);
				MC.clear();
				MC.reserve(key + 1 + current_clique_size);
				for(ui k = j;k < vs_size;k ++) MC.pb(vs[k]);
				for(ui k = current_clique_size;k > 0;k --) MC.pb(current_clique[k-1]);

#ifndef NDEBUG
				ui total_edges = MC.size();
				for(ui i = 0;i < MC.size();i ++) vis[MC[i]] = 1;
				for(ui i = 0;i < MC.size();i ++) for(ept j = pstart_o[MC[i]];j < pstart_o[MC[i]+1];j ++) {
					if(vis[edges_o[j]]) total_edges += 2;
				}
				if(total_edges != MC.size()*MC.size()) printf("WA! Not a clique\n");
				for(ui i = 0;i < MC.size();i ++) vis[MC[i]] = 0;
#endif
			}
			break;
		}
		vis[v] = 0;
		if(sparse) {
			for(ui k = pstart[v];k < pend[v];k ++) if(vis[edges[k]] == 1) heap->decrement(edges[k], 1);
		}
		else {
			for(ui k = pstart[v];k < pend[v];k ++) if(vis[edges[k]] == 1) -- degree[edges[k]];
		}
	}
	for(ui j = 0;j < vs_size;j ++) vis[vs[j]] = 0;

	return start_color;
}

//construct a better maximal clique and a tighter upper bound by considering all ego-networks
ui SubSearcher::ego_degen(const ui *peel_sequence, const ui *core, const ui *color, ui *local_UBs, const ui UB) {
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(!vis[i]);
	assert(MC.size() >= 2&&MC.size() <= n);
#endif

	// Timer t;

	ui max_local_UB = 0, initial_size = MC.size();
	ui *queue = new ui[maxCore];
	ListLinearHeap *heap = new ListLinearHeap(n, maxCore);
	heap->init(0, 0, nullptr, nullptr);
	for(ui i = n;i > 0;i --) {
		ui u = peel_sequence[i-1];
		local_UBs[i-1] = 0;

		if(n-i+1 <= MC.size()) continue;
		if(core[u] < MC.size()) break;

		//get N^+(u)
		ui vs_size = 0;
		get_higher_neighbors(u, vs_size, vs_buf, color_buf, pstart_o, edges_o);
		assert(vs_size <= maxCore);

		//color-based prune
		if(vs_size < MC.size()||color_bound(vs_buf.data(), vs_size, color, vis) < MC.size()) continue;

		for(ui j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 1;

		//test max_clique \cup \{u\}
		if(MC.size() > initial_size) greedy_extend(u, vis);

		//construct G[N^+(u)] and reduce by k-core
		ui old_size = vs_size, original_size = old_size;
		construct_induced_subgraph(vs_buf.data(), vs_size, vis, degree, pstart_o, edges_o);
		kcore_reduction(vs_buf.data(), vs_size, vis, degree, MC.size()-1, queue);
		for(ui j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 0;
		if(vs_size < old_size&&color_bound(vs_buf.data(), vs_size, color, vis) < MC.size()) continue;

		char kernel = 0; // by default set as 0
		ui current_clique_size = 1; current_clique[0] = u;
		if(kernel) {
			//construct matrix
			ui *rdegree = degree;
			construct_matrix(vs_buf.data(), vs_size, mapping, vis, rdegree, pstart_o, edges_o);

			contractions.clear();
			old_size = vs_size;
			degree_one_two_reduction_with_folding_matrix(current_clique_size, vs_buf.data(), vs_size, rdegree, contractions);
			if(MC.size() >= UB) break;

			//if(vs_size + current_clique_size <= max_clique.size()) continue;
			if(vs_size < old_size&&current_clique_size + color_bound(vs_buf.data(), vs_size, mapping, color, vis) <= MC.size()) continue;

			//degeneracy-based maximal clique
			//sort(vs_buf.begin(), vs_buf.begin()+vs_size);
			for(ui j = 0;j < vs_size;j ++) degree[vs_buf[j]] = vs_size - 1 - rdegree[vs_buf[j]];
			ui start_color = degeneracy_maximal_clique_matrix(current_clique_size, vs_buf.data(), vs_size, degree, 1, 0);
			if(MC.size() >= UB) break;

			//color-based upper bound
			//printf("start_color: %u, vs_size: %u\n", start_color, vs_size);
			local_UBs[i-1] = current_clique_size + coloring_matrix(vs_buf.data(), vs_size, rid, vis, 0, start_color);
		}
		else {
			//degeneracy-based maximal clique
			ui start_color = degeneracy_maximal_clique_adj_list(current_clique_size, vs_buf.data(), vs_size, vis, degree, heap);

			if(MC.size() >= UB) break;

			//color-based upper bound
			local_UBs[i-1] = current_clique_size + coloring_csr(vs_buf.data(), vs_size, original_size, rid, vis, 0, start_color);
		}

		if(local_UBs[i-1] > max_local_UB) max_local_UB = local_UBs[i-1];
	}

	delete[] queue;
	delete heap;

	if(max_local_UB > UB) max_local_UB = UB;

	ui new_UB = MC.size();
	if(max_local_UB > new_UB) new_UB = max_local_UB;

	// printf("*** ego_degen clique size: %lu, UB: %u, Time: %s (microseconds)\n", MC.size(), new_UB, Utility::integer_to_string(t.elapsed()).c_str());

	return new_UB;
}

void SubSearcher::heuristic_max_clique_max_degree(ui processed_threshold) {
	// Timer t;
	ui *head = new ui[n];
	ui *next = new ui[n];
	ui *degree = new ui[n];

	ui *vis = new ui[n];
	memset(vis, 0, sizeof(ui)*n);

	int max_degree = 0;
	for(ui i = 0;i < n;i ++) head[i] = n;
	for(ui i = 0;i < n;i ++) {
		degree[i] = pstart[i+1]-pstart[i];
		if(degree[i] > max_degree) max_degree = degree[i];
		next[i] = head[degree[i]];
		head[degree[i]] = i;
	}

	for(ui processed_vertices = 0;max_degree >= MC.size()&&processed_vertices < processed_threshold;processed_vertices ++) {
		ui u = n;
		while(max_degree >= MC.size()&&u == n) {
			for(ui v = head[max_degree];v != n;) {
				ui tmp = next[v];
				if(degree[v] == max_degree) {
					u = v;
					head[max_degree] = tmp;
					break;
				}
				else if(degree[v] >= MC.size()) {
					next[v] = head[degree[v]];
					head[degree[v]] = v;
				}
				v = tmp;
			}
			if(u == n) {
				head[max_degree] = n;
				-- max_degree;
			}
		}
		if(u == n) break;

		vis[u] = 1;
		for(ui k = pstart[u];k < pstart[u+1];k ++) if(!vis[edges[k]]) -- degree[edges[k]];

		vector<ui> vs;
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(!vis[edges[j]]) vs.pb(edges[j]);

		vector<ui> vs_deg(vs.size());
		for(ui j = 0;j < vs.size();j ++) vis[vs[j]] = 2;
		for(ui j = 0;j < vs.size();j ++) {
			ui v = vs[j], d = 0;
			for(ui k = pstart[v];k < pstart[v+1];k ++) {
				if(vis[edges[k]] == 2) ++ d;
			}
			vs_deg[j] = d;
		}
		for(ui j = 0;j < vs.size();j ++) vis[vs[j]] = 0;

		vector<ui> res; res.pb(u);
		ui vs_size = vs.size();
		while(vs_size > 0&&res.size() + vs_size > MC.size()) {
			ui idx = 0;
			for(ui j = 1;j < vs_size;j ++) {
				if(vs_deg[j] > vs_deg[idx]) idx = j;
				else if(vs_deg[j] == vs_deg[idx]&&degree[vs[j]] > degree[vs[idx]]) idx = j;
			}
			u = vs[idx];

			ui new_size = 0;
			for(ui j = pstart[u];j < pstart[u+1];j ++) if(!vis[edges[j]]) vis[edges[j]] = 2;
			for(ui j = 0;j < vs_size;j ++) if(vis[vs[j]]) {
				if(j != new_size) swap(vs[new_size], vs[j]);
				vs_deg[new_size] = vs_deg[j];
				++ new_size;
			}
			for(ui j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]] == 2) vis[edges[j]] = 0;

			res.pb(u);
			for(ui k = 0;k < new_size;k ++) vis[vs[k]] = k+2;
			for(ui j = new_size;j < vs_size;j ++) {
				ui v = vs[j];
				for(ui k = pstart[v];k < pstart[v+1];k ++) {
					if(vis[edges[k]] >= 2) -- vs_deg[vis[edges[k]]-2];
				}
			}
			for(ui k = 0;k < new_size;k ++) vis[vs[k]] = 0;

			vs_size = new_size;
		}

		if(res.size() > MC.size()) MC = res;
	}

	delete[] vis;
	delete[] head;
	delete[] next;
	delete[] degree;

	// printf("*** Heuristic clique size: %lu, time: %s (microseconds)\n", MC.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

//Build the oriented graph
void SubSearcher::shrink_graph(ui *&seq, ui *&core, ui *&color, ui *&out_mapping, ept *&pstart_o, ui *&edges_o) {
	ui *pos = new ui[n]; for(ui i = 0;i < n;i ++) pos[seq[i]] = i;
	ui cliqueSz = MC.size(); pstart_o = new ept[n+1]; pstart_o[0] = 0;
	for(ui i = 0;i < n;i ++) {
		pstart_o[i+1] = pstart_o[i];
		if(core[i] < cliqueSz) continue;
		for(ept j = pstart[i];j < pstart[i+1];j ++) if(pos[edges[j]] > pos[i]) edges[pstart_o[i+1] ++] = edges[j];
	}
	ui cnt = 0;
	for(ui i = 0;i < n;i ++) if(core[i] >= cliqueSz) {
		core[cnt] = core[i]; color[cnt] = color[i];
		out_mapping[cnt] = i; pos[i] = cnt ++;
	}
	ui x = n - cnt; for(ui i = x;i < n;i ++) seq[i-x] = pos[seq[i]];

	edges_o = new ui[pstart_o[n]];
	for(ui i = 0;i < pstart_o[n];i ++) edges_o[i] = pos[edges[i]];
	for(ui i = 0;i < cnt;i ++) pstart_o[i] = pstart_o[out_mapping[i]];
	pstart_o[cnt] = pstart_o[n];
	n = cnt;

	delete[] pos;

	ui *t_seq = new ui[n];
	memcpy(t_seq, seq, sizeof(ui)*n);
	delete[] seq; seq = t_seq;

	ui *t_core = new ui[n];
	memcpy(t_core, core, sizeof(ui)*n);
	delete[] core; core = t_core;

	ui *t_color = new ui[n];
	memcpy(t_color, color, sizeof(ui)*n);
	delete[] color; color = t_color;

	ept *t_pstart = new ept[n+1];
	memcpy(t_pstart, pstart_o, sizeof(ept)*(n+1));
	delete[] pstart_o; pstart_o = t_pstart;

	// printf("\tReduced graph size: |V|=%s, |E|=%s (undirected)\n", Utility::integer_to_string(cnt).c_str(), Utility::integer_to_string(pstart_o[n]).c_str());
}

void SubSearcher::recursive_search_clique_color_with_kernelization(const ui level, ui clique_size, ui vs_begin, ui vs_size, char kernel) {
	assert(clique_size <= MC.size()&&vs_buf.size() == color_buf.size()&&vs_size);
	ui *vs = vs_buf.data() + vs_begin;
	ui *color = color_buf.data() + vs_begin;

#ifndef _KERNEL_
	printf("Wrong invocation of recursive_search_clique_color_with_kernelization!\n");
#endif

#ifdef _STATISTIC_
	++ branches;
	if(level > max_depth) max_depth = level;
#endif

#ifndef NDEBUG
	if(!kernel) {
		ui color_cnt = 1;
		for(ui i = 1;i < vs_size;i ++) if(color[i] != color[i-1]) {
			++ color_cnt;
			for(ui j = i+1;j < vs_size;j ++) assert(color[j] != color[i-1]);
		}
		assert(color_cnt + clique_size == MC.size()+1);
	}
	for(ui i = 0;i < vs_size;i ++) for(ui j = 1;j < clique_size;j ++) assert(vs[i] != current_clique[j]);
	for(ui i = 0;i < vs_size;i ++) for(ui j = i+1;j < vs_size;j ++) assert(vs[i] != vs[j]);
	for(ui i = 1;i < clique_size;i ++) for(ui j = i+1;j < clique_size;j ++) assert(current_clique[i] == n || current_clique[i] != current_clique[j]);
#endif

	if(clique_size + 3 > MC.size()) {
		if(kernel) search_triangle_matrix(vs, vs_size, clique_size);
		else search_triangle_matrix_color(vs, vs_size, color, clique_size);
		if(clique_size > MC.size()) store_a_larger_clique(clique_size, "search_triangle", 1);
		return ;
	}

	kernel = 1;

	ui end_idx = 0, old_max_clique_size = MC.size();
	if(kernel) {
		obtain_degrees(vs, vs_size, degree);
		if(level) {
			ui *rdegree = degree;
			for(ui i = 0;i < vs_size;i ++) rdegree[vs[i]] = vs_size-1-degree[vs[i]];
			degree_one_two_three_reduction_with_folding_matrix(clique_size, vs, vs_size, rdegree, rid);
			for(ui i = 0;i < vs_size;i ++) degree[vs[i]] = vs_size-1-rdegree[vs[i]];
			if(clique_size > MC.size()) store_a_larger_clique(clique_size, "degree_one_two_three_with_folding", 1);

			if(!vs_size||MC.size() > old_max_clique_size) return ;
		}

		//ui cc_cnt = compute_connected_components(vs, vs_size);
		//if(cc_cnt > 1) ++ cc_larger_than_one[level];

		//degeneracy-based maximal clique
		ui start_color = degeneracy_maximal_clique_matrix(clique_size, vs, vs_size, degree, 1, 1);

		if(MC.size() > old_max_clique_size) return ;

		ui threshold = MC.size() - clique_size;
		//printf("start_color: %u, vs_size: %u\n", start_color, vs_size);
		ui color_cnt = coloring_matrix_advanced(vs, vs_size, color, start_color, threshold);
		if(color_cnt <= threshold) return ;

		//reduce to (threshold+1)-core
		if(reduce(vs, vs_size, color, threshold, clique_size)) return ;

		//reorganize vs
		end_idx = split_vs(vs, vs_size, color, threshold);

		if(color_cnt <= threshold+1) kernel = 0;
	}
	else {
		move_min_cardinality_color_to_front(vs, vs_size, color);

		end_idx = 1;
		while(end_idx < vs_size&&color[end_idx] == color[0]) ++ end_idx;
	}

	while(vs_buf.size() < vs_begin+vs_size+vs_size) {
		vs_buf.pb(0);
		color_buf.pb(0);
	}
	vs = vs_buf.data() + vs_begin;
	color = color_buf.data() + vs_begin;

	for(ui i = end_idx;i > 0&&MC.size() == old_max_clique_size;i --) {
		ui *tvs = vs + vs_size;
		ui *tcolor = color + vs_size;
		ui tvs_end = 0;
		unsigned char *t_matrix = matrix + vs[i-1]*matrix_len;

		ui upper_bound = 0;
		for(ui j = (kernel?i:end_idx);j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
			assert(vs_buf.size() > vs_begin + vs_size + tvs_end);

			if(tvs_end == 0||color[j] != tcolor[tvs_end-1]) ++ upper_bound;

			tvs[tvs_end] = vs[j];
			tcolor[tvs_end ++] = color[j];
		}
		// Note that a vertex u may not connected to a vertex with color i for i < color(u)
		if(upper_bound + 1 + clique_size <= MC.size()) continue;

		current_clique[clique_size] = vs[i-1];
		ui new_clique_size = clique_size + 1;
		ui old_size = tvs_end;

		if(!kernel) {
			assert(upper_bound + clique_size == MC.size());
			if(kernelization_color(new_clique_size, tvs, tvs_end, tcolor)) continue;
		}

		ui old_contraction_size = contractions.size(), old_changes_size = changes.size();
		recursive_search_clique_color_with_kernelization(level+1, new_clique_size, vs_begin + vs_size, tvs_end, kernel);
		while(contractions.size() > old_contraction_size) contractions.pop_back();
		while(changes.size() > old_changes_size) {
			std::pair<ui, ui> p = changes.back(); changes.pop_back();
			reverse_bit(matrix+p.first*matrix_len, p.second);
			reverse_bit(matrix+p.second*matrix_len, p.first);
		}
		if(!kernel&&old_size + end_idx == vs_size) break;

		vs = vs_buf.data() + vs_begin;
		color = color_buf.data() + vs_begin;
	}
}

void SubSearcher::recursive_search_clique_color_without_kernelization(const ui level, ui clique_size, ui vs_begin, ui vs_size) {
	assert(clique_size <= MC.size()&&vs_buf.size() == color_buf.size()&&vs_size);
	ui *vs = vs_buf.data() + vs_begin;
	ui *color = color_buf.data() + vs_begin;

#ifdef _STATISTIC_
	++ branches;
	if(level > max_depth) max_depth = level;
#endif

#ifndef NDEBUG
	//printf("vs_begin: %u\n", vs_begin);
	for(ui i = 0;i < vs_size;i ++) for(ui j = 1;j < clique_size;j ++) assert(vs[i] != current_clique[j]);
	for(ui i = 0;i < vs_size;i ++) for(ui j = i+1;j < vs_size;j ++) assert(vs[i] != vs[j]);
	for(ui i = 1;i < clique_size;i ++) for(ui j = i+1;j < clique_size;j ++) assert(current_clique[i] == n||current_clique[i] != current_clique[j]);
#endif

	if(clique_size == MC.size()) {
		current_clique[clique_size ++] = vs[0];
		store_a_larger_clique(clique_size, "search", 0);
		return ;
	}

	ui old_max_clique_size = MC.size();
	obtain_degrees(vs, vs_size, degree);

	//degeneracy_ordering
	ui start_color = degeneracy_maximal_clique_matrix(clique_size, vs, vs_size, degree, 1, 0);
	//printf("start_color: %u, vs_size: %u\n", start_color, vs_size);
	if(MC.size() > old_max_clique_size) return ;


	//printf("start_color: %u, vs_size: %u\n", start_color, vs_size);

	ui threshold = MC.size() - clique_size;
	ui color_cnt = coloring_matrix_advanced(vs, vs_size, color, start_color, threshold);
	if(color_cnt <= threshold) return ;

	//reorganize vs
	ui end_idx = split_vs(vs, vs_size, color, threshold);

	while(vs_buf.size() < vs_begin+vs_size+vs_size) {
		vs_buf.pb(0);
		color_buf.pb(0);
	}
	vs = vs_buf.data() + vs_begin;
	color = color_buf.data() + vs_begin;

	for(ui i = end_idx;i > 0&&MC.size() == old_max_clique_size;i --) {
		ui *tvs = vs + vs_size;
		ui *tcolor = color + vs_size;
		ui tvs_end = 0;
		unsigned char *t_matrix = matrix + vs[i-1]*matrix_len;

		ui upper_bound = 0;
		for(ui j = i;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
			assert(vs_buf.size() > vs_begin + vs_size + tvs_end);

			if(tvs_end == 0||color[j] != tcolor[tvs_end-1]) ++ upper_bound;

			tvs[tvs_end] = vs[j];
			tcolor[tvs_end ++] = color[j];
		}
		// Note that a vertex u may not connected to a vertex with color i for i < color(u)
		if(upper_bound + 1 + clique_size <= MC.size()) continue;

		current_clique[clique_size] = vs[i-1];
		ui new_clique_size = clique_size + 1;

		recursive_search_clique_color_without_kernelization(level+1, new_clique_size, vs_begin + vs_size, tvs_end);

		vs = vs_buf.data() + vs_begin;
		color = color_buf.data() + vs_begin;
	}
}

void SubSearcher::search_oriented(const ui *seq, const ui *core, const ui *color, const ui *subUB) {
#ifdef _STATISTIC_
	ui matrix_cnt = 0;
	double total_density = 0;
	double min_density = 1;
	long total_kernel_effect = 0;

	ui search_ego_cnt = 0;
	double total_ego_density = 0;
	double min_ego_density = 1;

	branches = 0;
	max_depth = 0;
#endif

	ui initial_size = MC.size();
	ui *queue = new ui[maxCore];

	for(ui i = n;i > 0;i --) {
		ui u = seq[i-1];
		if(subUB[i-1] <= MC.size()) continue;
		if(core[u] < MC.size()) break;
		ui oldSize = MC.size();

		//get N^+(u)
		ui vs_size = 0;
		get_higher_neighbors(u, vs_size, vs_buf, color_buf, pstart_o, edges_o);

		//color-based prune
		if(vs_size < MC.size()||color_bound(vs_buf.data(), vs_size, color, vis) < MC.size()) continue;

		for(ui j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 1;

		//test max_clique \cup \{u\}
		if(MC.size() > initial_size && greedy_extend(u, vis)) {
			for(ui j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 0;
			continue;
		}

#ifdef _KERNEL_
		//construct G[N^+(u)] and reduce by k-core
		ui old_size = vs_size;
		construct_induced_subgraph(vs_buf.data(), vs_size, vis, degree, pstart_o, edges_o);
		kcore_reduction(vs_buf.data(), vs_size, vis, degree, MC.size()-1, queue);
		for(ui j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 0;
		if(vs_size < old_size&&color_bound(vs_buf.data(), vs_size, color, vis) < MC.size()) continue;
#endif

		//construct matrix
		ui *rdegree = degree;
		construct_matrix(vs_buf.data(), vs_size, mapping, vis, rdegree, pstart_o, edges_o);

#ifdef _STATISTIC_
		++ matrix_cnt;
		ui total_edges = 0;
		for(ui j = 0;j < vs_size;j ++) total_edges += vs_size - 1 - rdegree[vs_buf[j]];
		double density = double(total_edges)/(vs_size*(vs_size-1));
		if(density < min_density) min_density = density;
		total_density += density;
#endif

		ui current_clique_size = 1; current_clique[0] = u;
		contractions.clear(); changes.clear();

#ifdef _KERNEL_
		old_size = vs_size;
		//degree_one_two_reduction_with_folding_matrix(current_clique_size, vs_buf.data(), vs_size, rdegree, contractions);
		degree_one_two_three_reduction_with_folding_matrix(current_clique_size, vs_buf.data(), vs_size, rdegree, rid);
		if(current_clique_size > MC.size()) store_a_larger_clique(current_clique_size, "outside kernelization", 1);
		// total_kernel_effect += old_size - vs_size;
		if(MC.size() > oldSize||!vs_size) continue;

#ifdef _STATISTIC_
		++ search_ego_cnt;
		total_edges = 0;
		for(ui j = 0;j < vs_size;j ++) total_edges += vs_size - 1 - rdegree[vs_buf[j]];
		density = double(total_edges)/(vs_size*(vs_size-1));
		if(density < min_ego_density) min_ego_density = density;
		total_ego_density += density;
#endif

		changes.clear();
		recursive_search_clique_color_with_kernelization(0, current_clique_size, 0, vs_size, 1);
		//recursive_search_clique_color_without_kernelization(current_clique_size, 0, vs_size);
#else
		recursive_search_clique_color_without_kernelization(0, current_clique_size, 0, vs_size);
#endif
	}

	delete[] queue;

#ifdef _STATISTIC_
	if(matrix_cnt == 0) printf("No matrix is constructed!\n");
	else {
		printf("Number of matrix constructed: %s\n", Utility::integer_to_string(matrix_cnt).c_str());
		printf("Average density: %.4lf, min density: %.4lf, average kernel_effect: %.4lf\n", total_density/matrix_cnt, min_density, double(total_kernel_effect)/matrix_cnt);

		printf("Number of egos searched: %s, branches: %s\n", Utility::integer_to_string(search_ego_cnt).c_str(), Utility::integer_to_string(branches).c_str());
		if(search_ego_cnt == 0) search_ego_cnt = 1;
		printf("Average ego_density: %.4lf, min ego_density: %.4lf\n", total_ego_density/search_ego_cnt, min_ego_density);
	}
#endif
}

ui SubSearcher::color_bound(const ui *vs, const ui vs_size, const ui *color, char *vis) {
#ifndef NDEBUG
	for(ui i = 0;i < vs_size;i ++) assert(!vis[color[vs[i]]]);
#endif
	ui color_bound = 0;
	for(ui i = 0;i < vs_size;i ++) if(!vis[color[vs[i]]]) {
		vis[color[vs[i]]] = 1;
		++ color_bound;
	}
	for(ui i = 0;i < vs_size;i ++) vis[color[vs[i]]] = 0;
	return color_bound;
}

ui SubSearcher::color_bound(const ui *vs, const ui vs_size, const ui *mapping, const ui *color, char *vis) {
#ifndef NDEBUG
	for(ui i = 0;i < vs_size;i ++) assert(!vis[color[mapping[vs[i]]]]);
#endif
	ui color_bound = 0;
	for(ui i = 0;i < vs_size;i ++) if(!vis[color[mapping[vs[i]]]]) {
		vis[color[mapping[vs[i]]]] = 1;
		++ color_bound;
	}
	for(ui i = 0;i < vs_size;i ++) vis[color[mapping[vs[i]]]] = 0;
	return color_bound;
}

//coloring a graph that is represented by matrix
//return the number of colors used
ui SubSearcher::coloring_matrix(const ui *vs, const ui vs_size, ui *color, char *vis, const ui start_idx, const ui start_color) {
	assert(start_color > 0&&start_color <= vs_size);
	for(ui i = vs_size - start_color;i < vs_size;i ++) color[vs[i]] = vs_size - i - 1;
	ui max_color = start_color-1;

	for(ui i = vs_size - start_color;i > start_idx;i --) {
		ui u = vs[i-1];
		unsigned char *t_matrix = matrix + u*matrix_len;
		for(ui j = i;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) vis[color[vs[j]]] = 1;
		for(ui j = 0;;j ++) if(!vis[j]) {
			color[u] = j;
			if(j > max_color) max_color = j;
			break;
		}
		for(ui j = i;j < vs_size;j ++) vis[color[vs[j]]] = 0;
	}

	return max_color + 1;
}

//coloring a graph that is represented by matrix, aiming to minimize the number of vertices with color >= threshold
//return the number of colors used
ui SubSearcher::coloring_matrix_advanced(const ui *vs, const ui vs_size, ui *color, const ui start_color, const ui threshold) {
	assert(rid != nullptr&&head != nullptr); // rid is used to temporarily store the color
	for(ui i = 0;i < vs_size;i ++) head[i] = n;

	assert(start_color > 0&&start_color <= vs_size);
	for(ui i = vs_size - start_color;i < vs_size;i ++) {
		ui c = vs_size - i - 1;
		rid[vs[i]] = c;
		id[i] = vs[i];
		next[i] = head[c];
		head[c] = i;
	}

	for(ui i = vs_size - start_color;i > 0;) {
		ui u = vs[-- i];
		unsigned char *t_matrix = matrix + u*matrix_len;
		rid[u] = n;
		for(ui j = 0;j < threshold;j ++) {
			char ok = 1;
			for(ui k = head[j];k != n;k = next[k]) if(test_bit(t_matrix, id[k])) {
				ok = 0;
				break;
			}
			if(ok) {
				rid[u] = j;
				id[i] = vs[i];
				next[i] = head[j];
				head[j] = i;
				break;
			}
		}
#ifndef _RECOLOR_
		continue;
#endif
		if(rid[u] < threshold) continue;

		for(ui j = 0;j < threshold;j ++) {
			ui cnt = 0, idx;
			for(ui k = head[j];k != n;k = next[k]) if(test_bit(t_matrix, id[k])) {
				++ cnt;
				if(cnt == 1) idx = k;
				else break;
			}
			assert(cnt > 0);
			if(cnt != 1) continue;
			unsigned char *tt_matrix = matrix + id[idx]*matrix_len;
			for(ui ii = threshold;ii > 0;) {
				-- ii;
				if(ii == j) continue;

				char ok = 1;
				for(ui k = head[ii];k != n;k = next[k]) if(test_bit(tt_matrix, id[k])) {
					ok = 0;
					break;
				}
				if(ok) {
					rid[id[idx]] = ii;
					id[i] = id[idx];
					next[i] = head[ii];
					head[ii] = i;

					rid[u] = j;
					id[idx] = vs[i];
					break;
				}
			}
			if(rid[u] < threshold) break;
		}
		if(rid[u] < threshold) continue;

		for(ui j = threshold;;j ++) {
			char ok = 1;
			for(ui k = head[j];k != n;k = next[k]) if(test_bit(t_matrix, id[k])) {
				ok = 0;
				break;
			}
			if(ok) {
				rid[u] = j;
				id[i] = vs[i];
				next[i] = head[j];
				head[j] = i;
				break;
			}
		}
	}

	ui max_color = 0;
	for(ui i = 0;i < vs_size;i ++) {
		color[i] = rid[vs[i]];
		if(color[i] > max_color) max_color = color[i];
	}

#ifndef NDEBUG
	for(ui i = 0;i < vs_size;i ++) {
		unsigned char *t_matrix = matrix + vs[i]*matrix_len;
		for(ui j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) assert(color[i] != color[j]);
	}
#endif
	return max_color+1;
}

// construct induced subgraph from pstart_o and edges_o
void SubSearcher::construct_induced_subgraph(const ui *vs, const ui vs_size, char *vis, ui *degree, const ept *pstart_o, const ui *edges_o) {
#ifndef NDEBUG
	for(ui i = 0;i < vs_size;i ++) assert(vis[vs[i]]);
#endif
	for(ui j = 0;j < vs_size;j ++) degree[vs[j]] = 0;
	for(ui j = 0;j < vs_size;j ++) {
		ui v = vs[j];
		for(ui k = pstart_o[v];k < pstart_o[v+1];k ++) {
			if(vis[edges_o[k]]) {
				++ degree[v];
				++ degree[edges_o[k]];
			}
		}
	}

	pstart[vs[0]] = 0;
	for(ui j = 1;j < vs_size;j ++) pstart[vs[j]] = pstart[vs[j-1]] + degree[vs[j-1]];
	for(ui j = 0;j < vs_size;j ++) pend[vs[j]] = pstart[vs[j]];
	for(ui j = 0;j < vs_size;j ++) {
		ui v = vs[j];
		for(ui k = pstart_o[v];k < pstart_o[v+1];k ++) {
			ui w = edges_o[k];
			if(vis[w]) {
				edges[pend[v] ++] = w;
				edges[pend[w] ++] = v;
			}
		}
	}
}

// construct matrix from pstart_o and edges_o for vertices in vs
void SubSearcher::construct_matrix(ui *vs, const ui vs_size, ui *mapping, char *vis, ui *rdegree, const ept *pstart_o, const ui *edges_o) {
	assert(rid != nullptr&&mapping != nullptr&&matrix != nullptr);
	assert(vs_size <= maxCore);

	mapping_n = vs_size;
#ifdef _BITSET_
	matrix_len = (vs_size+7)/8;
#else
	matrix_len = vs_size;
#endif
	for(ui j = 0;j < vs_size;j ++) {
		vis[vs[j]] = 1;
		rid[vs[j]] = j;
		mapping[j] = vs[j];
		vs[j] = j;
	}

	memset(matrix, 0, sizeof(unsigned char)*mapping_n*matrix_len);
	for(ui j = 0;j < mapping_n;j ++) rdegree[j] = mapping_n - 1;
	for(ui j = 0;j < mapping_n;j ++) for(ui k = pstart_o[mapping[j]];k < pstart_o[mapping[j]+1];k ++) if(vis[edges_o[k]]) {
		ui v = rid[edges_o[k]];
		assert(v >= 0&&v < mapping_n);
		set_bit(matrix + j*matrix_len, v);
		set_bit(matrix + v*matrix_len, j);
		-- rdegree[j]; -- rdegree[v];
	}
	for(ui j = 0;j < mapping_n;j ++) set_bit(matrix + j*matrix_len, j);

	for(ui j = 0;j < mapping_n;j ++) vis[mapping[j]] = 0;
}

ui SubSearcher::degeneracy_maximal_clique_matrix(ui current_clique_size, ui *vs, const ui vs_size, ui *degree, char heuristic_gen, char print) {
#ifndef NDEBUG
	for(ui i = 0;i < vs_size;i ++) {
		ui d = 0;
		unsigned char *t_matrix = matrix + vs[i]*matrix_len;
		for(ui j = 0;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) ++ d;
		assert(test_bit(t_matrix, vs[i]));
		assert(degree[vs[i]]+1 == d);
	}
#endif

	ui start_color = 0;
	for(ui j = 0;j < vs_size;j ++) {
		ui min_idx = j;
		for(ui k = j+1;k < vs_size;k ++) if(degree[vs[k]] < degree[vs[min_idx]]) min_idx = k;
		if(min_idx != j) std::swap(vs[min_idx], vs[j]);
		ui v = vs[j];
		if(degree[v] + 1 + j == vs_size) {
			for(ui k = j+1;k < vs_size;k ++) degree[vs[k]] = vs_size-1-k;
			start_color = degree[v];
			assert(start_color > 0&&start_color < vs_size);
			if(heuristic_gen&&degree[v] + 1 + current_clique_size > MC.size()) {
				for(ui k = j;k < vs_size;k ++) current_clique[current_clique_size ++] = vs[k];
				store_a_larger_clique(current_clique_size, "degen", print);
			}
			break;
		}
		unsigned char *t_matrix = matrix + v*matrix_len;
		for(ui k = j+1;k < vs_size;k ++) if(test_bit(t_matrix, vs[k])) -- degree[vs[k]];
	}
	return start_color;
}

void SubSearcher::degree_one_two_reduction_with_folding_matrix(ui &current_clique_size, ui *vs, ui &vs_size, ui *rdegree, std::vector<ui> &contractions) {
#ifndef NDEBUG
	for(ui i = 0;i < vs_size;i ++) {
		ui rd = 0;
		unsigned char *t_matrix = matrix + vs[i]*matrix_len;
		for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) ++ rd;
		assert(test_bit(t_matrix, vs[i]));
		assert(rdegree[vs[i]] == rd);
	}
#endif

	char changed = 1;
	while(changed) {
		changed = 0;
		for(ui i = 0;i < vs_size;) {
			ui rd = rdegree[vs[i]];
			if(rd == 0) {
				changed = 1;
				current_clique[current_clique_size ++] = vs[i];
				vs[i] = vs[-- vs_size];
			}
			else if(vs_size - rd + current_clique_size <= MC.size()) {
				changed = 1;
				remove_vertex_idx(vs, i, vs_size, rdegree);
			}
			else if(rd == 1) {
				changed = 1;
				current_clique[current_clique_size ++] = vs[i];
				unsigned char *t_matrix = matrix + vs[i]*matrix_len;
				vs[i] = vs[-- vs_size];
				for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
					remove_vertex_idx(vs, j, vs_size, rdegree);
					break;
				}
			}
			else if(rd == 2) {
				changed = 1;
				ui u = vs[i];
				unsigned char *t_matrix = matrix + u*matrix_len;
				vs[i] = vs[-- vs_size];
				ui idx1 = n, idx2 = n;
				for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
					if(idx1 == n) idx1 = j;
					else {
						idx2 = j;
						break;
					}
				}
				assert(idx1 != n&&idx2 != n&&idx1 < idx2);

				if(!test_bit(matrix+vs[idx1]*matrix_len, vs[idx2])) {
					current_clique[current_clique_size ++] = u;
					remove_vertex_idx(vs, idx2, vs_size, rdegree);
					remove_vertex_idx(vs, idx1, vs_size, rdegree);
				}
				else {
					current_clique[current_clique_size ++] = n;
					-- rdegree[vs[idx1]];
					contractions.pb(vs[idx1]); contractions.pb(vs[idx2]); contractions.pb(u); contractions.pb(1);
					unsigned char *t_matrix1 = matrix + vs[idx1]*matrix_len;
					unsigned char *t_matrix2 = matrix + vs[idx2]*matrix_len;
					vs[idx2] = vs[-- vs_size];

					for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix2, vs[j])) {
						if(test_bit(t_matrix1, vs[j])) {
							reverse_bit(t_matrix1, vs[j]);
							reverse_bit(matrix + vs[j]*matrix_len, vs[idx1]);
							++ rdegree[vs[idx1]];
						}
						else -- rdegree[vs[j]];
					}
				}
			}
			else ++ i;
		}
	}

#ifndef NDEBUG
	for(ui i = 0;i < vs_size;i ++) {
		unsigned char *t_matrix = matrix + vs[i]*matrix_len;
		ui rd = 0;
		for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) ++ rd;
		assert(rd == rdegree[vs[i]]);
		assert(vs_size - 1 - rd + current_clique_size >= MC.size());
		assert(rd >= 2);
	}
#endif

	if(current_clique_size > MC.size()) store_a_larger_clique(current_clique_size, "degree_one_two_with_folding", 0);
}

void SubSearcher::put_into_one_vector(ui &current_clique_size, ui &i, ui *vs, ui &vs_size, bool check, const ui rd, ui *rid) {
	if(vs_size - rd + current_clique_size <= MC.size()) {
		del.pb(vs[i]);
		//printf("insert %u\n", vs[i]);
	}
	else if(check&&rd <= 3) {
		switch(rd) {
		case 0: current_clique[current_clique_size ++] = vs[i];
				vs[i] = vs[-- vs_size]; rid[vs[i]] = i;
				-- i;
				break;
		case 1: degree_one.pb(vs[i]);
				break;
		case 2: degree_two.pb(vs[i]);
				break;
		case 3: degree_three.pb(vs[i]);
		}
	}
}

void SubSearcher::put_into_one_vector_eq(ui &current_clique_size, ui &i, ui *vs, ui &vs_size, bool check, const ui rd, ui *rid) {
	if(check) {
		if(rd <= 3&&vs_size - rd + current_clique_size > MC.size()) {
			switch(rd) {
			case 0: current_clique[current_clique_size ++] = vs[i];
					//printf("added %u to maximum clique\n", vs[i]);
					vs[i] = vs[-- vs_size]; rid[vs[i]] = i;
					-- i;
					break;
			case 1: degree_one.pb(vs[i]);
				break;
			case 2: degree_two.pb(vs[i]);
				break;
			case 3: degree_three.pb(vs[i]);
			}
		}
	}
	else {
		if(vs_size - rd + current_clique_size == MC.size()) {
			del.pb(vs[i]);
			//printf("insert %u\n", vs[i]);
		}
	}
}

void SubSearcher::degree_one_two_three_reduction_with_folding_matrix(ui &current_clique_size, ui *vs, ui &vs_size, ui *rdegree, ui *rid) {
	assert(del.empty()&&degree_one.empty()&&degree_two.empty()&&degree_three.empty());

	for(ui i = 0;i < vs_size;i ++) {
		rid[vs[i]] = i;
		put_into_one_vector(current_clique_size, i, vs, vs_size, true, rdegree[vs[i]], rid);
	}

	while(!del.empty()||!degree_one.empty()||!degree_two.empty()||!degree_three.empty()) {
		while(!del.empty()) {
			ui u = del.back(); del.pop_back();
			rdegree[u] = 0;
			//printf("u: %u, rid[u]: %u, vs[rid[u]]: %u, vs_size: %u\n", u, rid[u], vs[rid[u]], vs_size);
			assert(rid[u] < vs_size&&vs[rid[u]] == u);
			ui idx = rid[u];
			vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;
			unsigned char *t_matrix = matrix + u*matrix_len;
			for(ui i = 0;i < vs_size;i ++) {
				ui &rd = rdegree[vs[i]];
				ui old_rd = rd;
				if(!test_bit(t_matrix, vs[i])) -- rd;
				put_into_one_vector_eq(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
			}

#ifndef NDEBUG
			for(ui i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
		}
		while(del.empty()&&!degree_one.empty()) {
			ui u = degree_one.back(); degree_one.pop_back();
			if(rdegree[u] != 1) continue;

			current_clique[current_clique_size ++] = u;
			assert(rid[u] < vs_size&&vs[rid[u]] == u);
			ui idx = rid[u]; rdegree[u] = 0;
			vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;

			unsigned char *t_matrix = matrix + u*matrix_len; idx = n;
			for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
				idx = j;
				break;
			}
			assert(idx != n);

			t_matrix = matrix + vs[idx]*matrix_len;
			rdegree[vs[idx]] = 0;
			vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;
			for(ui i = 0;i < vs_size;i ++) {
				ui &rd = rdegree[vs[i]];
				ui old_rd = rd;
				if(!test_bit(t_matrix, vs[i])) -- rd;
				put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
			}

#ifndef NDEBUG
			for(ui i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
		}
		while(del.empty()&&degree_one.empty()&&!degree_two.empty()) {
			ui u = degree_two.back(); degree_two.pop_back();
			if(rdegree[u] != 2) continue;

			unsigned char *t_matrix = matrix + u*matrix_len;
			assert(rid[u] < vs_size&&vs[rid[u]] == u);
			ui idx = rid[u]; rdegree[u] = 0;
			vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;

			ui idx1 = n, idx2 = n;
			for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
				if(idx1 == n) idx1 = j;
				else {
					idx2 = j;
					break;
				}
			}
			assert(idx1 != n&&idx2 != n&&idx1 < idx2);

			unsigned char *t_matrix1 = matrix + vs[idx1]*matrix_len;
			unsigned char *t_matrix2 = matrix + vs[idx2]*matrix_len;
			if(!test_bit(t_matrix1, vs[idx2])) {
				current_clique[current_clique_size ++] = u;
				rdegree[vs[idx2]] = 0;
				vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
				rdegree[vs[idx1]] = 0;
				vs[idx1] = vs[-- vs_size]; rid[vs[idx1]] = idx1;

				for(ui i = 0;i < vs_size;i ++) {
					ui &rd = rdegree[vs[i]];
					ui old_rd = rd;
					if(!test_bit(t_matrix1, vs[i])) -- rd;
					if(!test_bit(t_matrix2, vs[i])) -- rd;

					put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
				}
			}
			else {
				current_clique[current_clique_size ++] = n;
				ui v = vs[idx1], old_degree_v = rdegree[v];
				-- rdegree[v];
				contractions.pb(v); contractions.pb(vs[idx2]); contractions.pb(u); contractions.pb(1); //type-1 contraction

				rdegree[vs[idx2]] = 0;
				vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;

				for(ui i = 0;i < vs_size;i ++) if(vs[i] != v) {
					ui &rd = rdegree[vs[i]];
					ui old_rd = rd;

					if(!test_bit(t_matrix2, vs[i])) {
						if(test_bit(t_matrix1, vs[i])) {
							reverse_bit(t_matrix1, vs[i]);
							reverse_bit(matrix + vs[i]*matrix_len, v);
							++ rdegree[v];
							changes.pb(std::make_pair(v, vs[i]));
						}
						else -- rd;
					}
					put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
				}
				put_into_one_vector(current_clique_size, rid[v], vs, vs_size, rdegree[v] != old_degree_v, rdegree[v], rid);
			}

#ifndef NDEBUG
			for(ui i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
		}

		while(del.empty()&&degree_one.empty()&&degree_two.empty()&&!degree_three.empty()) {
			ui u = degree_three.back(); degree_three.pop_back();
			if(rdegree[u] != 3) continue;

			unsigned char *t_matrix = matrix + u*matrix_len;
			assert(rid[u] < vs_size&&vs[rid[u]] == u);
			ui idx = rid[u];
			//std::swap(vs[idx], vs[vs_size-1]); rid[vs[idx]] = idx; rid[vs[vs_size-1]] = vs_size - 1;
			rdegree[u] = 0;
			vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;

			ui idx1 = n, idx2 = n, idx3 = n;
			for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
				if(idx1 == n) idx1 = j;
				else if(idx2 == n) idx2 = j;
				else {
					idx3 = j;
					break;
				}
			}
			assert(idx1 != n&&idx2 != n&&idx3 != n&&idx1 < idx2&&idx2 < idx3);

			unsigned char *t_matrix1 = matrix + vs[idx1]*matrix_len;
			unsigned char *t_matrix2 = matrix + vs[idx2]*matrix_len;
			unsigned char *t_matrix3 = matrix + vs[idx3]*matrix_len;
			char connected_1 = 0, connected_2 = 0, connected_3 = 0;
			if(test_bit(t_matrix1, vs[idx2])) connected_1 = 1;
			if(test_bit(t_matrix1, vs[idx3])) connected_2 = 1;
			if(test_bit(t_matrix2, vs[idx3])) connected_3 = 1;
			char total_connected = connected_1 + connected_2 + connected_3;
			if(total_connected == 0) { //isolation reduction
				//rdegree[u] = 0; -- vs_size;

				current_clique[current_clique_size ++] = u;
				rdegree[vs[idx3]] = 0;
				vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;
				rdegree[vs[idx2]] = 0;
				vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
				rdegree[vs[idx1]] = 0;
				vs[idx1] = vs[-- vs_size]; rid[vs[idx1]] = idx1;

				for(ui i = 0;i < vs_size;i ++) {
					ui &rd = rdegree[vs[i]];
					ui old_rd = rd;
					if(!test_bit(t_matrix1, vs[i])) -- rd;
					if(!test_bit(t_matrix2, vs[i])) -- rd;
					if(!test_bit(t_matrix3, vs[i])) -- rd;

					put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
				}
			}
			else if(total_connected == 1) {
				//rdegree[u] = 0; -- vs_size;

				current_clique[current_clique_size ++] = n;
				//the following ensures that vs[idx2] is the vertex to be deleted
				if(connected_3) {
					std::swap(vs[idx1], vs[idx2]);
					std::swap(t_matrix1, t_matrix2);
					rid[vs[idx1]] = idx1;
				}
				else if(connected_1) {
					std::swap(vs[idx2], vs[idx3]);
					std::swap(t_matrix2, t_matrix3);
				}

				ui v = vs[idx1], old_degree_v = rdegree[v];
				rdegree[v] -= 2;
				contractions.pb(v); contractions.pb(vs[idx3]); contractions.pb(u); contractions.pb(1); //type-1 contraction

				rdegree[vs[idx3]] = 0; rdegree[vs[idx2]] = 0;
				assert(idx2 != idx3);
				if(idx2 < idx3) {
					vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;
					vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
				}
				else {
					vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
					vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;
				}

				for(ui i = 0;i < vs_size;i ++) if(vs[i] != v) {
					ui &rd = rdegree[vs[i]];
					ui old_rd = rd;

					if(!test_bit(t_matrix2, vs[i])) -- rd;

					if(!test_bit(t_matrix3, vs[i])) {
						if(test_bit(t_matrix1, vs[i])) {
							reverse_bit(t_matrix1, vs[i]);
							reverse_bit(matrix + vs[i]*matrix_len, v);
							changes.pb(std::make_pair(v, vs[i]));
							++ rdegree[v];
						}
						else -- rd;
					}
					put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
				}
				put_into_one_vector(current_clique_size, rid[v], vs, vs_size, rdegree[v] != old_degree_v, rdegree[v], rid);
			}
			else if(total_connected == 2) {
				//rdegree[u] = 0; -- vs_size;

				current_clique[current_clique_size ++] = n;
				//the following ensures that vs[idx2] is the vertex to be deleted
				if(!connected_3) {
					std::swap(vs[idx1], vs[idx3]);
					std::swap(t_matrix1, t_matrix3);
					rid[vs[idx1]] = idx1;
				}
				else if(!connected_2) {
					std::swap(vs[idx2], vs[idx3]);
					std::swap(t_matrix2, t_matrix3);
					rid[vs[idx2]] = idx2;
				}

				ui v = vs[idx1], old_degree_v = rdegree[v];
				ui w = vs[idx2], old_degree_w = rdegree[w];
				-- rdegree[v]; -- rdegree[w];
				contractions.pb(v); contractions.pb(w); contractions.pb(vs[idx3]); contractions.pb(u); contractions.pb(2); //type-2 contraction

				rdegree[vs[idx3]] = 0;
				vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;

				for(ui i = 0;i < vs_size;i ++) if(vs[i] != v&&vs[i] != w) {
					ui &rd = rdegree[vs[i]];
					ui old_rd = rd;

					if(!test_bit(t_matrix3, vs[i])) {
						-- rd;
						if(test_bit(t_matrix1, vs[i])) {
							reverse_bit(t_matrix1, vs[i]);
							reverse_bit(matrix + vs[i]*matrix_len, v);
							changes.pb(std::make_pair(v, vs[i]));
							++ rdegree[v];
							++ rd;
						}

						if(test_bit(t_matrix2, vs[i])) {
							reverse_bit(t_matrix2, vs[i]);
							reverse_bit(matrix + vs[i]*matrix_len, w);
							changes.pb(std::make_pair(w, vs[i]));
							++ rdegree[w];
							++ rd;
						}
					}
					put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
				}
				put_into_one_vector(current_clique_size, rid[v], vs, vs_size, rdegree[v] != old_degree_v, rdegree[v], rid);
				put_into_one_vector(current_clique_size, rid[w], vs, vs_size, rdegree[w] != old_degree_w, rdegree[w], rid);
			}
			else {
				assert(total_connected == 3);
				//rdegree[u] = 0; -- vs_size;

				current_clique[current_clique_size ++] = n;

				ui v1 = vs[idx1], old_degree_v1 = rdegree[v1];
				ui v2 = vs[idx2], old_degree_v2 = rdegree[v2];
				ui v3 = vs[idx3], old_degree_v3 = rdegree[v3];
				-- rdegree[v3];
				reverse_bit(t_matrix1, v2); reverse_bit(t_matrix2, v1);
				changes.pb(std::make_pair(v1, v2));
				contractions.pb(v1); contractions.pb(v2); contractions.pb(v3); contractions.pb(u); contractions.pb(3); //type-3 contraction

				for(ui i = 0;i < vs_size;i ++) if(vs[i] != v1&&vs[i] != v2&&vs[i] != v3) {
					ui &rd = rdegree[vs[i]];
					ui old_rd = rd;

					char conn1 = test_bit(t_matrix1, vs[i]);
					char conn2 = test_bit(t_matrix2, vs[i]);
					char conn3 = test_bit(t_matrix3, vs[i]);

					if(conn1&&!conn2) {
						reverse_bit(t_matrix1, vs[i]);
						reverse_bit(matrix+vs[i]*matrix_len, v1);
						changes.pb(std::make_pair(v1, vs[i]));
						++ rd;
						++ rdegree[v1];
					}

					if(conn2&&!conn3) {
						reverse_bit(t_matrix2, vs[i]);
						reverse_bit(matrix+vs[i]*matrix_len, v2);
						changes.pb(std::make_pair(v2, vs[i]));
						++ rd;
						++ rdegree[v2];
					}

					if(conn3&&!conn1) {
						reverse_bit(t_matrix3, vs[i]);
						reverse_bit(matrix+vs[i]*matrix_len, v3);
						changes.pb(std::make_pair(v3, vs[i]));
						++ rd;
						++ rdegree[v3];
					}

					put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
				}
				put_into_one_vector(current_clique_size, rid[v1], vs, vs_size, rdegree[v1] != old_degree_v1, rdegree[v1], rid);
				put_into_one_vector(current_clique_size, rid[v2], vs, vs_size, rdegree[v2] != old_degree_v2, rdegree[v2], rid);
				put_into_one_vector(current_clique_size, rid[v3], vs, vs_size, rdegree[v3] != old_degree_v3, rdegree[v3], rid);
			}

#ifndef NDEBUG
			for(ui i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
		}
	}

#ifndef NDEBUG
	for(ui i = 0;i < vs_size;i ++) assert(rdegree[vs[i]] > 3&&vs_size - rdegree[vs[i]] + current_clique_size > MC.size());
#endif
}

void SubSearcher::get_higher_neighbors(const ui u, ui &vs_size, std::vector<ui> &vs_buf, std::vector<ui> &color_buf, const ept *pstart_o, const ui *edges_o) {
	vs_size = 0;
	for(ui j = pstart_o[u];j < pstart_o[u+1];j ++) {
		if(vs_buf.size() == vs_size) {
			vs_buf.pb(edges_o[j]);
			color_buf.pb(0);
		}
		else vs_buf[vs_size] = edges_o[j];
		++ vs_size;
	}
}

//greedily enlarge max_clique by including u if feasible
char SubSearcher::greedy_extend(const ui u, const char *vis) {
	for(auto v: MC) if(!vis[v]) return 0; MC.pb(u);
	return 1;
}

void SubSearcher::kcore_reduction(ui *vs, ui &vs_size, char *vis, ui *degree, const ui K, ui *queue) {
#ifndef NDEBUG
	for(ui i = 0;i < vs_size;i ++) assert(vis[vs[i]]);
#endif

	ui queue_n = 0;
	for(ui i = 0;i < vs_size;i ++) if(degree[vs[i]] < K) {
		vis[vs[i]] = 0;
		queue[queue_n ++] = vs[i];
	}
	for(ui i = 0;i < queue_n;i ++) {
		ui u = queue[i];
		for(ui j = pstart[u];j < pend[u];j ++) if(vis[edges[j]]) {
			if((-- degree[edges[j]]) == K-1) {
				vis[edges[j]] = 0;
				queue[queue_n ++] = edges[j];
			}
		}
	}
	ui new_vs_size = 0;
	for(ui i = 0;i < vs_size;i ++) if(vis[vs[i]]) {
		if(new_vs_size != i) std::swap(vs[new_vs_size], vs[i]);
		++ new_vs_size;
	}
	vs_size = new_vs_size;
}

char SubSearcher::kernelization_color(ui &clique_size, ui *vs, ui &vs_size, ui *color) {
#ifndef NDEBUG
	for(ui i = 0;i < vs_size-1;i ++) if(color[i] != color[i+1]) {
		for(ui j = i+1;j < vs_size;j ++) assert(color[i] != color[j]);
	}
#endif
	while(true) {
		ui idx = n;
		for(ui i = 0;i < vs_size;i ++) {
			if((i==0||color[i]!=color[i-1])&&(i+1==vs_size||color[i]!=color[i+1])) {
				idx = i;
				break;
			}
		}
		if(idx == n) break;

		current_clique[clique_size ++] = vs[idx];

		ui new_size = 0;
		unsigned char *t_matrix = matrix + vs[idx]*matrix_len;
		for(ui i = 0;i < vs_size;i ++) {
			if(i == idx) continue;
			if(test_bit(t_matrix, vs[i])) {
				vs[new_size] = vs[i];
				color[new_size ++] = color[i];
			}
			else if(i+1 == vs_size||color[i+1] != color[i]) {
				if(new_size == 0||color[i] != color[new_size-1]) return 1;
			}
		}
		vs_size = new_size;
	}

	if(clique_size > MC.size()) {
		store_a_larger_clique(clique_size, "kernel_color", 1);
		return 1;
	}
	return 0;
}

void SubSearcher::move_min_cardinality_color_to_front(ui *vs, const ui vs_size, ui *color) {
	ui c, min_cnt = n, cnt = 0;
	for(ui i = 0;i < vs_size;i ++) {
		++ cnt;
		if(i+1 == vs_size||color[i] != color[i+1]) {
			if(cnt < min_cnt) {
				min_cnt = cnt;
				c = color[i];
			}
			cnt = 0;
		}
	}
	ui idx = vs_size;
	while(idx > 0&&color[idx-1] != c) -- idx;
	for(ui i = idx;i > 0;i --) if(color[i-1] != c) {
		std::swap(vs[i-1], vs[-- idx]);
		std::swap(color[i-1], color[idx]);
	}
}

void SubSearcher::obtain_degrees(const ui *vs, const ui vs_size, ui *degree) {
	for(ui i = 0;i < vs_size;i ++) degree[vs[i]] = 0;
	for(ui i = 0;i < vs_size;i ++) {
		ui &d = degree[vs[i]];
		unsigned char *t_matrix = matrix + vs[i]*matrix_len;
		for(ui j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
			++ degree[vs[j]];
			++ d;
		}
	}
}

char SubSearcher::reduce(ui *vs, ui &vs_size, ui *color, const ui threshold, ui clique_size) {
#ifndef NDEBUG
	for(ui i = 0;i < vs_size;i ++) {
		ui d = 0;
		unsigned char *t_matrix = matrix + vs[i]*matrix_len;
		for(ui j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) ++ d;
		assert(degree[vs[i]] == d);
	}
#endif

	assert(rid != nullptr);
	ui start = 0;
	for(;start < vs_size;start ++) if(degree[vs[start]] >= threshold) {
		ui *neighbors = rid;
		ui neighbors_n = 0;
		unsigned char *t_matrix = matrix + vs[start]*matrix_len;
		ui color_cnt = 0;
		for(ui i = start+1;i < vs_size;i ++) if(test_bit(t_matrix, vs[i])) {
			neighbors[neighbors_n ++] = i;
			if(!vis[color[i]]) {
				vis[color[i]] = 1;
				++ color_cnt;
			}
		}
		assert(neighbors_n == degree[vs[start]]);
		for(ui i = 0;i < neighbors_n;i ++) {
			vis[color[neighbors[i]]] = 0;
			neighbors[i] = vs[neighbors[i]];
		}
		if(color_cnt < threshold) continue;

		if(degree[vs[start]] >= threshold + 1) break;

		char ok = 1;
		for(ui i = 0;i < neighbors_n&&ok;i ++) {
			t_matrix = matrix + neighbors[i]*matrix_len;
			for(ui j = i+1;j < neighbors_n;j ++) if(!test_bit(t_matrix, neighbors[j])) {
				ok = 0;
				break;
			}
		}
		if(ok) {
			current_clique[clique_size ++] = vs[start];
			for(ui i = 0;i < neighbors_n;i ++) current_clique[clique_size ++] = neighbors[i];
			assert(clique_size == MC.size()+1);
			store_a_larger_clique(clique_size, "reduce", 1);
			return 1;
		}
	}

	if(start) {
		for(ui i = start;i < vs_size;i ++) {
			vs[i-start] = vs[i];
			color[i-start] = color[i];
		}
		vs_size -= start;
	}

	if(vs_size == 0) return 1;
	return 0;
}

void SubSearcher::remove_vertex_idx(ui *vs, ui idx, ui &vs_size, ui *rdegree) {
	assert(idx >= 0&&idx < vs_size);

	unsigned char *t_matrix = matrix + vs[idx]*matrix_len;
	vs[idx] = vs[-- vs_size];
	for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) -- rdegree[vs[j]];
}

void SubSearcher::search_triangle_matrix(const ui *vs, const ui vs_size, ui &clique_size) {
	if(clique_size+2 == MC.size()) {
		for(ui i = 0;i < vs_size&&clique_size < MC.size();i ++) {
			unsigned char *t_matrix1 = matrix + vs[i]*matrix_len;
			assert(rid != nullptr);
			ui *neighbors = rid;
			ui neighbors_n = 0;
			for(ui j = i+1;j < vs_size;j ++) if(test_bit(t_matrix1, vs[j])) {
				neighbors[neighbors_n ++] = vs[j];
			}
			for(ui j = 0;j < neighbors_n&&clique_size < MC.size();j ++) {
				unsigned char *t_matrix2 = matrix + neighbors[j]*matrix_len;
				for(ui k = j+1;k < neighbors_n;k ++) if(test_bit(t_matrix2, neighbors[k])) {
					current_clique[clique_size ++] = vs[i];
					current_clique[clique_size ++] = neighbors[j];
					current_clique[clique_size ++] = neighbors[k];
					break;
				}
			}
		}
	}
	else if(clique_size+1 == MC.size()) {
		for(ui i = 0;i < vs_size&&clique_size < MC.size();i ++) {
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(ui j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
				current_clique[clique_size ++] = vs[i];
				current_clique[clique_size ++] = vs[j];
				break;
			}
		}
	}
	else if(clique_size == MC.size()) current_clique[clique_size ++] = vs[0];
}

void SubSearcher::search_triangle_matrix_color(const ui *vs, const ui vs_size, const ui *color, ui &clique_size) {
	if(clique_size+2 == MC.size()) {
		ui idx1 = 1;
		while(color[idx1] == color[0]) ++ idx1;
		assert(idx1+1 < vs_size);
		ui idx2 = idx1+1;
		while(color[idx2] == color[idx1]) ++ idx2;
		assert(idx2 < vs_size);
		for(ui i = 0;i < idx1&&clique_size < MC.size();i ++) {
			unsigned char *t_matrix1 = matrix + vs[i]*matrix_len;
			assert(rid != nullptr);
			ui *neighbors = rid;
			ui neighbors_n = 0;
			for(ui j = idx2;j < vs_size;j ++) if(test_bit(t_matrix1, vs[j])) {
				neighbors[neighbors_n ++] = vs[j];
			}
			for(ui j = idx1;j < idx2&&clique_size < MC.size();j ++) if(test_bit(t_matrix1, vs[j])) {
				unsigned char *t_matrix2 = matrix + vs[j]*matrix_len;
				for(ui k = 0;k < neighbors_n;k ++) if(test_bit(t_matrix2, neighbors[k])) {
					current_clique[clique_size ++] = vs[i];
					current_clique[clique_size ++] = vs[j];
					current_clique[clique_size ++] = neighbors[k];
					break;
				}
			}
		}
	}
	else if(clique_size+1 == MC.size()) {
		ui idx1 = 1;
		while(color[idx1] == color[0]) ++ idx1;
		assert(idx1 < vs_size);
		for(ui i = 0;i < idx1&&clique_size < MC.size();i ++) {
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(ui j = idx1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
				current_clique[clique_size ++] = vs[i];
				current_clique[clique_size ++] = vs[j];
				break;
			}
		}
	}
	else if(clique_size == MC.size()) current_clique[clique_size ++] = vs[0];
}

ui SubSearcher::split_vs(ui *vs, const ui vs_size, ui *color, const ui threshold) {
	for(ui i = 0;i < threshold;i ++) head[i] = n;
	for(ui i = vs_size;i > 0;i --) if(color[i-1] < threshold) {
		ui j = i-1;
		next[vs[j]] = head[color[j]];
		head[color[j]] = vs[j];
	}
	ui end_idx = 0;
	for(ui i = 0;i < vs_size;i ++) if(color[i] >= threshold) {
		color[end_idx] = color[i];
		vs[end_idx] = vs[i];
		++ end_idx;
	}
	ui new_size = end_idx;
	for(ui i = threshold;i > 0;i --) for(ui u = head[i-1];u != n;u = next[u]) {
		vs[new_size] = u;
		color[new_size ++] = i-1;
	}
	assert(new_size == vs_size);

	return end_idx;
}

void SubSearcher::store_a_larger_clique(const ui clique_size, const char *info, char print) {
	assert(MC.size() < clique_size);
	while(MC.size() < clique_size) MC.pb(0);

	ui contract_size = contractions.size();
	for(ui k = clique_size-1;k > 0;k --) {
		if(current_clique[k] == n) {
			assert(contract_size > 0);
			if(contractions[contract_size-1] == 1) {
				contract_size -= 4;
				if(vis[contractions[contract_size]]) MC[k] = contractions[contract_size+1];
				else MC[k] = contractions[contract_size+2];
			}
			else if(contractions[contract_size-1] == 2) {
				contract_size -= 5;
				assert(!vis[contractions[contract_size]]||!vis[contractions[contract_size+1]]);
				if(!vis[contractions[contract_size]]&&!vis[contractions[contract_size+1]]) MC[k] = contractions[contract_size+3];
				else MC[k] = contractions[contract_size+2];
			}
			else if(contractions[contract_size-1] == 3) {
				contract_size -= 5;
				char in1 = (vis[contractions[contract_size]] != 0);
				char in2 = (vis[contractions[contract_size+1]] != 0);
				char in3 = (vis[contractions[contract_size+2]] != 0);
				assert(in1+in2+in3 < 3);
				if(in1+in2+in3 == 0) MC[k] = contractions[contract_size+3];
				else if(in1+in2+in3 == 1) {
					if(in1) MC[k] = contractions[contract_size+1];
					else if(in2) MC[k] = contractions[contract_size+2];
					else MC[k] = contractions[contract_size];
				}
				else {
					if(!in1) MC[k] = contractions[contract_size];
					else if(!in2) MC[k] = contractions[contract_size+1];
					else MC[k] = contractions[contract_size+2];
				}
			}
			else {
				printf("WA!\n");
			}
		}
		else MC[k] = current_clique[k];
		vis[MC[k]] = 1;
	}
	assert(contract_size == 0);
	for(ui k = 1;k < MC.size();k ++) vis[MC[k]] = 0;
	for(ui k = 1;k < MC.size();k ++) MC[k] = mapping[MC[k]];
	MC[0] = current_clique[0];
	// if(print) printf("%s finds clique of size: %lu\n", info, MC.size());

#ifndef NDEBUG
	ui total_edges = MC.size();
	for(ui i = 0;i < MC.size();i ++) vis[MC[i]] = 1;
	for(ui i = 0;i < MC.size();i ++) for(ept j = pstart_o[MC[i]];j < pstart_o[MC[i]+1];j ++) {
		if(vis[edges_o[j]]) total_edges += 2;
	}
	if(total_edges != MC.size()*MC.size()) printf("WA! Not a clique! %s\n", info);
	for(ui i = 0;i < MC.size();i ++) vis[MC[i]] = 0;
#endif
}

void SubSearcher::release(){
	// CNei.clear();
	// CNeiRemap.clear();
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(edges != nullptr) {
		delete[] edges;
		edges = nullptr;
	}
	if(pend != nullptr) {
		delete[] pend;
		pend = nullptr;
	}
	if(pstart_o != nullptr) {
		delete[] pstart_o;
		pstart_o = nullptr;
	}
	if(edges_o != nullptr) {
		delete[] edges_o;
		edges_o = nullptr;
	}
	if(tmp_vs != nullptr) {
		delete[] tmp_vs;
		tmp_vs = nullptr;
	}
	if(tmp_color != nullptr) {
		delete[] tmp_color;
		tmp_color = nullptr;
	}
	if(matrix != nullptr) {
		delete[] matrix;
		matrix = nullptr;
	}
	if(head != nullptr) {
		delete[] head;
		head = nullptr;
	}
	if(next != nullptr) {
		delete[] next;
		next = nullptr;
	}
	if(vis != nullptr) {
		delete[] vis;
		vis = nullptr;
	}
	if(degree != nullptr) {
		delete[] degree;
		degree = nullptr;
	}
	if(rid != nullptr) {
		delete[] rid;
		rid = nullptr;
	}
	if(id != nullptr) {
		delete[] id;
		id = nullptr;
	}
	if(current_clique != nullptr) {
		delete[] current_clique;
		current_clique = nullptr;
	}
	if(mapping != nullptr) {
		delete[] mapping;
		mapping = nullptr;
	}
}
