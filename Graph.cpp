#include "Graph.h"
#include "MSearcher.h"

using namespace std;
Graph::Graph(const char *_dir, const int _K) {
	dir = string(_dir);
	K = _K;

	n = m = 0;

	pstart = nullptr;
	pend = pend_buf = nullptr;
	edges = nullptr;
	maxDeg=0;
	KDC.clear();

	s_degree = s_edges = NULL;
	s_pstart = s_pend = NULL;
	s_peel_sequence = s_core = NULL;
	s_vis = NULL;
	s_heap = NULL;

	s_edgelist_pointer = NULL;
	s_tri_cnt = s_edge_list = NULL;
	s_active_edgelist = NULL;
	s_deleted = NULL;
	
}

Graph::~Graph() {
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(pend != nullptr) {
		delete[] pend;
		pend = nullptr;
	}
	if(pend_buf != NULL) {
		delete[] pend_buf;
		pend_buf = NULL;
	}
	if(edges != nullptr) {
		delete[] edges;
		edges = nullptr;
	}
	if(s_degree != NULL) {
		delete[] s_degree;
		s_degree = NULL;
	}
	if(s_pstart != NULL) {
		delete[] s_pstart;
		s_pstart = NULL;
	}
	if(s_pend != NULL) {
		delete[] s_pend;
		s_pend = NULL;
	}
	if(s_edges != NULL) {
		delete[] s_edges;
		s_edges = NULL;
	}
	if(s_peel_sequence != NULL) {
		delete[] s_peel_sequence;
		s_peel_sequence = NULL;
	}
	if(s_core != NULL) {
		delete[] s_core;
		s_core = NULL;
	}
	if(s_vis != NULL) {
		delete[] s_vis;
		s_vis = NULL;
	}
	if(s_heap != NULL) {
		delete s_heap;
		s_heap = NULL;
	}
	if(s_edgelist_pointer != NULL) {
		delete[] s_edgelist_pointer;
		s_edgelist_pointer = NULL;
	}
	if(s_active_edgelist != NULL) {
		delete[] s_active_edgelist;
		s_active_edgelist = NULL;
	}
	if(s_deleted != NULL) {
		delete[] s_deleted;
		s_deleted = NULL;
	}
}

void Graph::read() {
	FILE *f = fopen(dir.c_str(), "rb");
	fread(&n, sizeof(int), 1, f); fread(&n, sizeof(int), 1, f); fread(&m, sizeof(int), 1, f);
#ifdef _TEST_
	printf("#n=%d\n#m=%d\n", n, m/2);
#else 
	printf("\tn = %d; m = %d (undirected)\n", n, m/2);
#endif
	

	int *degree = new int[n];
	fread(degree, sizeof(int), n, f);
	if(pstart == nullptr) pstart = new int[n+1];
	if(edges == nullptr) edges = new int[m];

	pstart[0] = 0;
	for(int i = 0;i < n;i ++) {
		if(degree[i] > 0) {
			
			fread(edges+pstart[i], sizeof(int), degree[i], f);

			// remove self loops and parallel edges
			int *buff = edges+pstart[i];
			sort(buff, buff+degree[i]);
			int idx = 0;
			for(int j = 0;j < degree[i];j ++) {
				if(buff[j] >= n) printf("vertex id %u wrong\n", buff[j]);
				if(buff[j] == i||(j > 0&&buff[j] == buff[j-1])) continue;//buff[j]==i-->exist selfloop, buff[j]==buff[j-1]-->exist parallel edges
				buff[idx ++] = buff[j];
			}
			degree[i] = idx;
			maxDeg=max(maxDeg, degree[i]);
		}

		pstart[i+1] = pstart[i] + degree[i];
	}
	printf("max degree: %d\n", maxDeg);
	fclose(f);
	delete[] degree;
}

void Graph::write() {
	FILE *fout = fopen("KDC.txt", "w");
	fprintf(fout, "%u\n", KDC.size());
	sort(KDC.begin(), KDC.end());
	for(int i = 0;i < KDC.size();i ++) fprintf(fout, "%d ", KDC[i]);
	fclose(fout);
}

void Graph::search() {
	Timer t;
	KDC.resize(K+1); //screen out trivial cases
	int *seq = new int[n];
	int *core = new int[n];
	int *deg = new int[n];
	char *vis = new char[n];

	ListLinearHeap *heap = new ListLinearHeap(n, n-1);
	int UB = degen(n, seq, core, pstart, edges, deg, vis, heap, true);
	delete heap;
	delete[] vis;
	delete[] deg;

	if(KDC.size() < UB) {		
		int old_size = KDC.size();
		int *out_mapping = new int[n];
		int *rid = new int[n];

		shrink_graph(n, m, seq, core, out_mapping, NULL, rid, pstart, edges);

		int *deg = new int[n]; for(int i = 0;i < n;i ++) deg[i] = pstart[i+1] - pstart[i];

		ListLinearHeap *linear_heap = new ListLinearHeap(n, n-1);
		linear_heap->init(n, n-1, seq, deg);

		pend = new int[n];
		pend_buf = new int[n];
		int *edge_list = new int[m];
		int *edgelist_pointer = new int[m];
		int *tri_cnt = new int[m/2];
		oriented_triangle_counting(n, m, seq, pstart, pend, edges, edgelist_pointer, rid); // edgelist_pointer currently stores triangle_counts
		reorganize_oriented_graph(n, tri_cnt, edge_list, pstart, pend, pend_buf, edges, edgelist_pointer, rid);
		int *active_edgelist = new int[m>>1]; int active_edgelist_n = m>>1;
		for(int i = 0;i < (m>>1);i ++) active_edgelist[i] = i;
		for(int i = 0;i < n;i ++) pend[i] = pstart[i+1];
		int *Qe = new int[m>>1];
		char *deleted = new char[m>>1];
		memset(deleted, 0, sizeof(char)*(m>>1));
		char *exists = new char[n];
		memset(exists, 0, sizeof(char)*n);
		int *Qv = new int[n]; int Qv_n = 0;
		
		m -= 2*peeling(n, linear_heap, Qv, Qv_n, KDC.size()-K, Qe, KDC.size()-K-1, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deg, pstart, pend, edges, exists);
		printf("*** Core-Truss Shrink: n = %d, m = %d, Density = %.2f%%\n", n-Qv_n, m/2, double(m)/(n-Qv_n)/(n-Qv_n-1)*100);

		Timer tt;

		int max_n = n - Qv_n;
		//This is the structure for subgraph
		s_degree = new int[max_n];
		s_pstart = new int[max_n+1];
		s_pend = new int[max_n];
		s_edges = new int[m];
		s_peel_sequence = new int[max_n];
		s_core = new int[max_n];
		s_vis = new char[max_n];
		s_heap = new ListLinearHeap(max_n,max_n-1);
		s_edgelist_pointer = new int[m];
		s_tri_cnt = new int[m/2];
		s_edge_list = new int[m];
		s_active_edgelist = new int[m/2];
		s_deleted = new char[m/2];

		vector<pair<int,int> > vp; vp.reserve(m/2);
		int *t_degree = new int[n];

		for(int i = 0;i < n&&m&&KDC.size() < UB;i ++) {
			int u, key; linear_heap->pop_min(u, key);
			if(key < KDC.size()-K) {
				if(deg[u] != 0) {
					Qv[0] = u; Qv_n = 1;
					m -= 2*peeling(n, linear_heap, Qv, Qv_n, KDC.size()-K, Qe, KDC.size()-K-1, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deg, pstart, pend, edges, exists);
				}
				continue;
			}
			if(!m) break;

			int *ids = Qv; int ids_n = 0; int subUB=UB;
			int pre_size0;
			do{
				pre_size0=KDC.size();
				induceSubgraph(u, ids, ids_n, rid, vp, Qe, t_degree, exists, pend, deleted, edgelist_pointer);
				if(ids_n > KDC.size()&& vp.size()*2 < m) subUB=min(UB,subgraph_heuri(ids, ids_n, vp, rid, Qv, Qe, exists));
			}
			while(KDC.size()!=pre_size0);
			// induceSubgraph(u, ids, ids_n, rid, vp, Qe, t_degree, exists, pend, deleted, edgelist_pointer);
			// if(ids_n > KDC.size()&& vp.size()*2 < m) subUB=min(UB,subgraph_heuri(ids, ids_n, vp, rid, Qv, Qe, exists));
			int pre_size = KDC.size();
			if(ids_n > pre_size && subUB> pre_size) { 
				if(subNum==4){
					int a=1;
				}
				maxSubSz=max(maxSubSz, ids_n);
				ExactSearcher exact_solver(K, ids_n, vp, KDC, UB); exact_solver.search();
				if(KDC.size() != pre_size) for(int j = 0;j < KDC.size();j ++) KDC[j] = ids[KDC[j]];
			}
			Qv[0] = u; Qv_n = 1;
			m -= 2*peeling(n, linear_heap, Qv, Qv_n, KDC.size()-K, Qe, KDC.size()-K-1, tri_cnt, active_edgelist, active_edgelist_n, edge_list, edgelist_pointer, deleted, deg, pstart, pend, edges, exists);
		}
		if(KDC.size() > old_size) for(int i = 0;i < KDC.size();i ++) KDC[i] = out_mapping[KDC[i]];
		
#ifdef _TEST_
		printf("#MaxKDCSize=%d\n#SearchTime=%.2f\n#TotalTime=%.2f\n", KDC.size(), double(tt.elapsed())/1000000, double(t.elapsed())/1000000);
		printf("#maxP=%d\n#minPUB=%d\n#maxME=%d\n", max_P_end, P_UBMin,maxME);
		printf("#NodeCount=%lld\n",tree_cnt);
#else
		printf("\tMaxKDC Size: %d, Search Time: %.2f, Total Time: %.2f\n", KDC.size(), double(tt.elapsed())/1000000, double(t.elapsed())/1000000);
		printf("max P_end: %d, min P UB: %d,max missing edges: %d, max Col Sz: %d\n", max_P_end, P_UBMin,maxME , maxColSz);
		printf("max sub-graph size: %d, avg subgraph dense: %lf\n",maxSubSz ,avgSubDense);
		printf("Max Mutual Exclusive: %d, Mutual exclusive sum: %lld, Mutual exclusive avg: %lf, Mutual exclusive dense avg: %lf\n",MaxMuExNum, MuExSum, (double)MuExSum/(tree_cnt+0.1), denSum/(denNum+0.001));
		printf("Sum of Mutual Exclusive in color set: %lld, max mutual exclusive in color set: %d\n", colMuExSum, maxColMuNum);
		printf("Search Tree Size: %lld, the sum of induced graphs: %lld, the num of sub graphs: %d, avg size: %.2f\n",tree_cnt, szSum, subNum, double(szSum)/subNum);
#endif
		
		
		if(tree_cnt>0){
#ifdef _TEST_
			printf("#Bound-Prune=%lld %.2f%%\n#Color-Bound-Prune=%lld %.2f%%\n", boundPrune,double(boundPrune)/(boundPrune+colorBndPrune+tree_cnt)*100, colorBndPrune, double(colorBndPrune)/(boundPrune+colorBndPrune+tree_cnt)*100);
#else
			printf("Bound Prune: %lld %.2f%%, Color Bound Prune: %lld %.2f%%\n", boundPrune,double(boundPrune)/(boundPrune+colorBndPrune+tree_cnt)*100, colorBndPrune, double(colorBndPrune)/(boundPrune+colorBndPrune+tree_cnt)*100);
#endif
		}
		delete linear_heap;
		delete[] t_degree;
		delete[] exists;
		delete[] out_mapping;
		delete[] rid;
		delete[] deg;
		delete[] edgelist_pointer;
		delete[] tri_cnt;
		delete[] active_edgelist;
		delete[] Qe;
		delete[] Qv;
		delete[] deleted;
	}else if(KDC.size()>=UB || UB<K+2){
#ifdef _TEST_
		printf("#MaxKDCSize=%d\n#SearchTime=%.2f\n#TotalTime=%.2f\n", KDC.size(), 0.0, double(t.elapsed())/1000000);
		printf("#NodeCount=%lld\n",tree_cnt);
#else
		printf("\tMaxKDC Size: %d, Search Time: %.2f, Total Time: %.2f\n", KDC.size(), 0.0, double(t.elapsed())/1000000);
		printf("Search Tree Size: %lld, the sum of induced graphs: %lld\n",0,0);
#endif
	}
	delete[] core;
	delete[] seq;

}

void Graph::load_graph_from_edgelist(int _n, const vector<pair<int,int> > &edge_list, int &n, int &m, int *degree, int *pstart, int *edges) {
	n = _n;
	m = (int)edge_list.size()*2;
	for(int i = 0; i < n; i++) degree[i] = 0;
	for(int i = 0;i < m/2;i ++) {
		assert(edge_list[i].first >= 0&&edge_list[i].first < n&&edge_list[i].second >= 0&&edge_list[i].second < n);
		degree[edge_list[i].first] ++;
		degree[edge_list[i].second] ++;
	}

	pstart[0] = 0;
	for(int i = 0;i < n;i ++) pstart[i+1] = pstart[i]+degree[i];
	for(int i = 0;i < m/2;i ++) {
		int a = edge_list[i].first, b = edge_list[i].second;
		edges[pstart[a]++] = b;
		edges[pstart[b]++] = a;
	}
	for(int i = 0;i < n;i ++) pstart[i] -= degree[i];
}

void Graph::extract_graph(int n, int m, int *degree, int *ids, int &ids_n, int *rid, vector<pair<int,int> > &vp, char *exists, int *pstart, int *pend, int *edges, char *deleted, int *edgelist_pointer) {
	ids_n = 0; vp.clear();
	for(int i=0; i<n; ++i){
		if(degree[i]){
			ids[ids_n] = i; rid[i] = ids_n++;
		}
	}
	for(int i = 0; i<ids_n ; i++) {
		int u = ids[i];
		for(int j = pstart[u];j < pend[u]; j++) if(!deleted[edgelist_pointer[j]] && u < edges[j]) {
			vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
	}
}

void Graph::extract_subgraph(int u, int *ids, int &ids_n, int *rid, vector<pair<int,int> > &vp, char *exists, int *pstart, int *pend, int *edges, char *deleted, int *edgelist_pointer) {
	ids_n = 0; vp.clear();
	ids[ids_n++] = u; exists[u] = 1; rid[u] = 0;
	int u_n = pstart[u];
	for(int i = pstart[u];i < pend[u];i ++) if(!deleted[edgelist_pointer[i]]) {
		edges[u_n] = edges[i]; edgelist_pointer[u_n++] = edgelist_pointer[i];
		int v = edges[i];
		rid[v] = ids_n; ids[ids_n++] = v; exists[v] = 1;
	}
	pend[u] = u_n;
	int old_size = ids_n;
	for(int i = 1;i < old_size;i ++) {
		u = ids[i];
		u_n = pstart[u];
		for(int j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			int v = edges[j];
			if(exists[v]) continue;
			rid[v] = ids_n; ids[ids_n++] = v; exists[v] = 1;
		}
		pend[u] = u_n;
	}
	for(int i = 0;i < old_size;i ++) {
		u = ids[i];
		for(int j = pstart[u];j < pend[u];j ++) if(edges[j] > u) {
			vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
	}
	for(int i = old_size;i < ids_n;i ++) {
		u = ids[i];
		u_n = pstart[u];
		for(int j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			if(edges[j] > u&&exists[edges[j]]) vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
		pend[u] = u_n;
	}
	for(int i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
}

void Graph::induceSubgraph(int u, int *ids, int &ids_n, int *rid, vector<pair<int,int> > &vp, int *Q, int* degree, char *exists, int *pend, char *deleted, int *edgelist_pointer) {
	vp.clear();
	ids_n = 0; ids[ids_n++] = u; exists[u] = 1;
	int u_n = pstart[u];
	for(int i = pstart[u];i < pend[u];i ++) if(!deleted[edgelist_pointer[i]]) {
		edges[u_n] = edges[i]; edgelist_pointer[u_n++] = edgelist_pointer[i];
		int v = edges[i];
		ids[ids_n++] = v; exists[v] = 2;
	}
	pend[u] = u_n;
	
	int Q_n = 0;
	for(int i = 1;i < ids_n;i ++) {
		u = ids[i];
		u_n = pstart[u];
		degree[u] = 0;
		for(int j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			if(exists[edges[j]] == 2) ++ degree[u];
		}
		pend[u] = u_n;
		if(degree[u]+2+K <= KDC.size()) Q[Q_n++] = u;
	}
	for(int i = 0;i < Q_n;i ++) {
		u = Q[i];
		exists[u] = 10;
		for(int j = pstart[u];j < pend[u];j ++) if(exists[edges[j]] == 2) {
			if( (--degree[edges[j]])+2+K == KDC.size()) {
				assert(Q_n < m/2);
				Q[Q_n++] = edges[j];
			}
		}
	}
	assert(Q_n <= ids_n);
	if(ids_n-Q_n+K<=KDC.size()) {
		for(int i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
		ids_n = 0;
		return ;
	}
	
	int nr_size = ids_n;
	for(int i = 1;i < nr_size;i ++) if(exists[ids[i]] == 2) {
		u = ids[i];
		for(int j = pstart[u];j < pend[u];j ++) {
			if(!exists[edges[j]]) {
				ids[ids_n++] = edges[j];
				exists[edges[j]] = 3;
				degree[edges[j]] = 1;
			}
			else if(exists[edges[j]] == 3) ++ degree[edges[j]];
		}
	}

#ifndef NDEBUG
	//printf("Entire list: ");
	//for(ui i = 0;i < nr_size;i ++) printf(" %u", ids[i]);
	//printf("\n");
#endif

	int new_size = 1;
	for(int i = 1;i < nr_size;i ++) {
		if(exists[ids[i]] == 10) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
#ifndef NDEBUG
	if(new_size + Q_n != nr_size) {
		printf("new_size: %u, Q_n: %u, nr_size: %u\n", new_size, Q_n, nr_size);
		printf("New list: ");
		for(int i = 0;i < new_size;i ++) printf(" %u", ids[i]);
		printf("\n");
		printf("Pruned list: ");
		for(int i = 0;i < Q_n;i ++) printf(" %u", Q[i]);
		printf("\n");
	}
#endif
	assert(new_size + Q_n == nr_size);
	int old_nr_size = nr_size;
	nr_size = new_size;
	for(int i = old_nr_size;i < ids_n;i ++) {
		if(degree[ids[i]]+2+K-1<= KDC.size()) exists[ids[i]] = 0;
		else ids[new_size++] = ids[i];
	}
	ids_n = new_size;
#ifndef NDEBUG
	assert(exists[ids[0]] == 1);
	for(int i = 1;i < nr_size;i ++) assert(exists[ids[i]] == 2);
	for(int i = nr_size;i < ids_n;i ++) assert(exists[ids[i]] == 3);
#endif

	//for(ui i = 0;i < ids_n;i ++) printf(" %u", ids[i]);
	//printf("\n");

	for(int i = 0;i < ids_n;i ++) {
		assert(exists[ids[i]]);
		rid[ids[i]] = i;
	}

	for(int i = 0;i < nr_size;i ++) {
		u = ids[i];
		for(int j = pstart[u];j < pend[u];j ++) if(exists[edges[j]]&&edges[j] > u) {
			assert(!deleted[edgelist_pointer[j]]);
			vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
	}
	for(int i = nr_size;i < ids_n;i ++) {
		u = ids[i];
		u_n = pstart[u];
		for(int j = pstart[u];j < pend[u];j ++) if(!deleted[edgelist_pointer[j]]) {
			edges[u_n] = edges[j]; edgelist_pointer[u_n++] = edgelist_pointer[j];
			if(edges[j] > u&&exists[edges[j]]) vp.push_back(make_pair(rid[u], rid[edges[j]]));
		}
		pend[u] = u_n;
	}
	for(int i = 0;i < ids_n;i ++) exists[ids[i]] = 0;
#ifndef NDEBUG
	for(int i = 0;i < n;i ++) assert(exists[i] == 0);
#endif
}

// degeneracy-based k-plex
// return an upper bound of the maximum k-plex size
int Graph::degen(int n, int *seq, int *core, int *pstart, int *edges, int *degree, char *vis, ListLinearHeap *heap, bool output) {
	Timer t;
	int threshold = KDC.size()-K;
	int edgeCnt=0;
	for(int i = 0;i < n;i++) degree[i] = pstart[i+1] - pstart[i], edgeCnt+=degree[i]; edgeCnt=edgeCnt/2;
	int queue_n = 0, new_size = 0;
	for(int i = 0;i < n;i++) if(degree[i] < threshold) seq[queue_n ++] = i;//first order reduction
	for(int i = 0;i < queue_n;i ++) {
		int u = seq[i]; degree[u] = 0;
		for(int j = pstart[u];j < pstart[u+1];j ++) if(degree[edges[j]] > 0) {
			if((degree[edges[j]] --) == threshold) seq[queue_n ++] = edges[j];
			edgeCnt--;
		}
	}
	int UB = n; if(queue_n == n) UB = KDC.size();
	memset(vis, 0, sizeof(char)*n);
	for(int i = 0;i < n;i ++) {
		if(degree[i] >= threshold) seq[queue_n + (new_size ++)] = i;
		else vis[i] = 1, core[i] = 0;
	}
	assert(queue_n + new_size == n);

	if(new_size != 0) {
		heap->init(new_size, new_size-1, seq+queue_n, degree);
		int maxCore = 0; int idx = n;
		int t_UB=0;
		for(int i = 0;i < new_size;i ++) {
			if(idx == n && edgeCnt + K >= (1LL*(new_size - i)*(new_size - i - 1)/2) ) idx = i;
			int u, key; heap->pop_min(u, key);
			if(key > maxCore) maxCore = key; core[u] = maxCore;
			seq[queue_n + i] = u; vis[u] = 1;
			t_UB=max(t_UB,min(core[u]+K+1,new_size-i));
			for(int j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]] == 0) {
				heap->decrement(edges[j], 1);
			}
			edgeCnt-=key;
		}
		UB=min(UB,min(t_UB,maxCore+(int)floor((1+sqrt(1+8*K))/2)));
		if(new_size - idx > KDC.size()) { 
			KDC.clear(); for(int i = idx;i < new_size;i ++) KDC.pb(seq[queue_n + i]);
			if(!output) printf("Find a KDC of size: %u\n", new_size - idx);
		}
#ifdef _TEST_
		if(output) printf("#HeuSize=%u\n#MaxCore=%d\n#UB=%d\n#HeuTime=%.2f\n", KDC.size(), maxCore, UB, double(t.elapsed())/1000000);
#else
		if(output) printf("*** HeuriKDC size: %u, MaxCore: %d, UB: %d, Heuri Time: %.2f\n", KDC.size(), maxCore, UB, double(t.elapsed())/1000000);
#endif
	}
	return UB;
}

// in_mapping and out_mapping can be the same array
// note that core is not maintained, and is assumed to not be used anymore
void Graph::shrink_graph(int &n, int &m, int *peel_sequence, int *core, int *out_mapping, int *in_mapping, int *rid, int *pstart, int *edges) {
	int cnt = 0;
	for(int i = 0;i < n;i ++) if(core[i] + K + 1> KDC.size()) {
		rid[i] = cnt;
		if(in_mapping == NULL) out_mapping[cnt] = i;
		else out_mapping[cnt] = in_mapping[i];
		++ cnt;
	}

	if(cnt != n) {
		cnt = 0;
		int pos = 0;
		for(int i = 0;i < n;i ++) if(core[i] + K + 1 > KDC.size()) {
			int t_start = pstart[i]; pstart[cnt] = pos;
			for(int j = t_start;j < pstart[i+1];j ++) if(core[edges[j]] + K + 1> KDC.size()) {
				edges[pos ++] = rid[edges[j]];
			}
			++ cnt;
		}
		pstart[cnt] = pos;

		//printf("%u %u %u %u\n", n, cnt, core[peel_sequence[n-cnt-1]], core[peel_sequence[n-cnt]]);
		assert(core[peel_sequence[n-cnt-1]] == 0||core[peel_sequence[n-cnt-1]] + K + 1<= KDC.size());
		assert(cnt == 0||core[peel_sequence[n-cnt]] + K + 1> KDC.size());
		for(int i = 0;i < cnt;i ++) {
			peel_sequence[i] = rid[peel_sequence[n-cnt+i]];
			//core[i] = core[out_mapping[i]];
		}

		n = cnt;
		m = pos;
	}

	printf("*** Core Shrink: n = %d, m = %d \n", n, m/2);
}

// orient graph and triangle counting
void Graph::oriented_triangle_counting(int n, int m, int *peel_sequence, int *pstart, int *pend, int *edges, int *tri_cnt, int *adj) {
	int *rid = adj;
	for(int i = 0;i < n;i ++) rid[peel_sequence[i]] = i;
	for(int i = 0;i < n;i ++) {
		int &end = pend[i] = pstart[i];
		for(int j = pstart[i];j < pstart[i+1];j ++) if(rid[edges[j]] > rid[i]) edges[end ++] = edges[j];
	}

#ifndef NDEBUG
	long long sum = 0;
	for(int i = 0;i < n;i ++) sum += pend[i] - pstart[i];
	// printf("%lld %lld\n", sum, m);
	assert(sum*2 == m);
#endif

	memset(adj, 0, sizeof(int)*n);
	long long cnt = 0;
	memset(tri_cnt, 0, sizeof(int)*m);
	for(int u = 0;u < n;u ++) {
		for(int j = pstart[u];j < pend[u];j ++) adj[edges[j]] = j+1;

		for(int j = pstart[u];j < pend[u];j ++) {
			int v = edges[j];
			for(int k = pstart[v];k < pend[v];k ++) if(adj[edges[k]]) {
				++ tri_cnt[j];
				++ tri_cnt[k];
				++ tri_cnt[adj[edges[k]]-1];
				++ cnt;
			}
		}

		for(int j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
	}
}

// reorganize the adjacency lists
// and sort each adjacency list to be in increasing order
void Graph::reorganize_oriented_graph(int n, int *tri_cnt, int *edge_list, int *pstart, int *pend, int *pend2, int *edges, int *edgelist_pointer, int *buf) {
	for(int i = 0;i < n;i ++) pend2[i] = pend[i];
	int pos = 0;
	for(int i = 0;i < n;i ++) {
		for(int j = pstart[i];j < pend[i];j ++) {
			tri_cnt[pos>>1] = edgelist_pointer[j]; edge_list[pos++] = i; edge_list[pos++] = edges[j];

			int &k = pend2[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j] = (pos>>1)-1;
			edges[k ++] = i;
		}
	}

#ifndef NDEBUG
	for(int i = 0;i < n;i ++) assert(pend2[i] == pstart[i+1]);
#endif

	for(int i = 0;i < n;i ++) {
		pend2[i] = pend[i];
		pend[i] = pstart[i];
	}
	for(int i = 0;i < n;i ++) {
		for(int j = pend2[i];j < pstart[i+1];j ++) {
			int &k = pend[edges[j]];
			edgelist_pointer[k] = edgelist_pointer[j];
			edges[k ++] = i;
		}
	}

	int *ids = pend2;
	for(int i = 0;i < n;i ++) {
		if(pend[i] == pstart[i]||pend[i] == pstart[i+1]) continue;
		int j = pstart[i], k = pend[i], pos = 0;
		while(j < pend[i]&&k < pstart[i+1]) {
			if(edges[j] < edges[k]) {
				ids[pos] = edges[j];
				buf[pos ++] = edgelist_pointer[j ++];
			}
			else {
				ids[pos] = edges[k];
				buf[pos ++] = edgelist_pointer[k ++];
			}
		}
		while(j < pend[i]) {
			ids[pos] = edges[j];
			buf[pos ++] = edgelist_pointer[j ++];
		}
		while(k < pstart[i+1]) {
			ids[pos] = edges[k];
			buf[pos ++] = edgelist_pointer[k ++];
		}
		for(int j = 0;j < pos;j ++) {
			edges[pstart[i] + j] = ids[j];
			edgelist_pointer[pstart[i] + j] = buf[j];
		}
	}
}


char Graph::find(int u, int w, int &b, int e, char *deleted, int &idx, int *edgelist_pointer, int *edges) {
	if(b >= e) return 0;

	while(b+1 < e) {
		idx = b + (e-b)/2;
		if(edges[idx] > w) e = idx;
		else b = idx;
	}

	if(edges[b] == w) {
		idx = edgelist_pointer[b];
		if(!deleted[idx]) return 1;
	}

	return 0;
}

// return the number of peeled edges
int Graph::peeling(int critical_vertex, ListLinearHeap *linear_heap, int *Qv, int &Qv_n, int d_threshold, int *Qe, int t_threshold, int *tri_cnt, int *active_edgelist, int &active_edgelist_n, int *edge_list, int *edgelist_pointer, char *deleted, int *degree, int *pstart, int *pend, int *edges, char *exists) {
	int Qe_n = 0;
#ifndef NO_TRUSS_PRUNE
	bool initialize_Qe=t_threshold>0;
	if(initialize_Qe) {
		int active_edgelist_newn = 0;
		for(int j = 0;j < active_edgelist_n;j ++) if(!deleted[active_edgelist[j]]) {
			if(tri_cnt[active_edgelist[j]] < t_threshold) Qe[Qe_n++] = active_edgelist[j];
			else active_edgelist[active_edgelist_newn ++] = active_edgelist[j];
		}
		active_edgelist_n = active_edgelist_newn;
	}
#endif

	//printf("%lu\n", Qe_n);

	int deleted_edges_n = 0;
	int Qv_idx = 0;
	while(Qv_idx < Qv_n || Qe_n) {
		if(Qe_n == 0) {
			//printf("hit\n");
			int u = Qv[Qv_idx ++]; // delete u from the graph due to have a degree < d_threshold
			int u_n = pstart[u];
			for(int k = pstart[u];k < pend[u];k ++) if(!deleted[edgelist_pointer[k]]) {
				edges[u_n] = edges[k]; edgelist_pointer[u_n++] = edgelist_pointer[k];
				exists[edges[k]] = 1;
			}
			pend[u] = u_n;

			for(int k = pstart[u];k < pend[u];k ++) deleted[edgelist_pointer[k]] = 1;
			deleted_edges_n += pend[u] - pstart[u];
			degree[u] = 0;
			if(linear_heap != NULL) linear_heap->del(u);
			//printf("Removed %u\n", u);

			for(int k= pstart[u];k < pend[u];k ++) {
				int v = edges[k];
				int v_n = pstart[v];
				for(int x = pstart[v];x < pend[v];x ++) if(!deleted[edgelist_pointer[x]]) {
					edges[v_n] = edges[x]; edgelist_pointer[v_n++] = edgelist_pointer[x];
					if(edges[x] > v&&exists[edges[x]]) {
						if( (tri_cnt[edgelist_pointer[x]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[x];
					}
				}
				pend[v] = v_n;

				if( (degree[v]--) == d_threshold) {
					Qv[Qv_n++] = v;
					if(v == critical_vertex) {
						for(int k = pstart[u];k < pend[u];k ++) exists[edges[k]] = 0;
						return 0;
					}
				}
				if(linear_heap != NULL) linear_heap->decrement(v, 1);
			}

			for(int k = pstart[u];k < pend[u];k ++) exists[edges[k]] = 0;
		}
#ifdef NO_TRUSS_PRUNE
		Qe_n = 0;
#endif
		for(int j = 0;j < Qe_n;j ++) {
			int idx = Qe[j];
			int u = edge_list[idx<<1], v = edge_list[(idx<<1)+1];
			int tri_n = tri_cnt[idx];
			//printf("remove %u %u\n", u, v);
			deleted[idx] = 1;
			if( (degree[u] --) == d_threshold) {
				Qv[Qv_n++] = u;
				if(u == critical_vertex) return 0;
			}
			if( (degree[v] --) == d_threshold) {
				Qv[Qv_n++] = v;
				if(v == critical_vertex) return 0;
			}
			//printf("before\n");
			if(linear_heap != NULL) {
				linear_heap->decrement(u, 1);
				linear_heap->decrement(v, 1);
			}
			//printf("after\n");
			deleted_edges_n ++;
			
			if(degree[u] < degree[v]) swap(u,v);
			//printf("here\n");

			if(degree[u] > degree[v]*2) { // binary search
			//if(false) {
				int v_n = pstart[v], start = pstart[u];
				for(int k = pstart[v];k < pend[v];k ++) if(!deleted[edgelist_pointer[k]]) {
					edges[v_n] = edges[k]; edgelist_pointer[v_n++] = edgelist_pointer[k];

					if(tri_n&&find(u, edges[k], start, pend[u], deleted, idx, edgelist_pointer, edges)) {
						-- tri_n;
						if( (tri_cnt[idx]--) == t_threshold) Qe[Qe_n++] = idx;
						if( (tri_cnt[edgelist_pointer[k]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[k];
					}
				}
				pend[v] = v_n;
				assert(tri_n == 0);
			}
			else { // sorted_merge
				int ii = pstart[u], jj = pstart[v];
				int u_n = pstart[u], v_n = pstart[v];

				while(true) {
					while(ii < pend[u]&&deleted[edgelist_pointer[ii]]) ++ ii;
					while(jj < pend[v]&&deleted[edgelist_pointer[jj]]) ++ jj;
					if(ii >= pend[u]||jj >= pend[v]) break;

					if(edges[ii] == edges[jj]) {
						edges[u_n] = edges[ii]; edgelist_pointer[u_n++] = edgelist_pointer[ii];
						edges[v_n] = edges[jj]; edgelist_pointer[v_n++] = edgelist_pointer[jj];

						if( (tri_cnt[edgelist_pointer[ii]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[ii];
						if( (tri_cnt[edgelist_pointer[jj]]--) == t_threshold) Qe[Qe_n++] = edgelist_pointer[jj];

						++ ii;
						++ jj;
					}
					else if(edges[ii] < edges[jj]) {
						edges[u_n] = edges[ii]; edgelist_pointer[u_n++] = edgelist_pointer[ii];
						++ ii;
					}
					else {
						edges[v_n] = edges[jj]; edgelist_pointer[v_n++] = edgelist_pointer[jj];
						++ jj;
					}
				}
				while(ii < pend[u]) {
					if(!deleted[edgelist_pointer[ii]]) {
						edges[u_n] = edges[ii];
						edgelist_pointer[u_n++] = edgelist_pointer[ii];
					}
					++ ii;
				}
				while(jj < pend[v]) {
					if(!deleted[edgelist_pointer[jj]]) {
						edges[v_n] = edges[jj];
						edgelist_pointer[v_n++] = edgelist_pointer[jj];
					}
					++ jj;
				}
				pend[u] = u_n; pend[v] = v_n;
			}
		}
		Qe_n = 0;
	}
	return deleted_edges_n;
}


int Graph::subgraph_heuri(ui *ids, ui &_n, vector<pair<int,int> > &edge_list, ui *rid, ui *Qv, ui *Qe, char *exists) {
	ui s_n; ept s_m;
	load_graph_from_edgelist(_n, edge_list, s_n, s_m, s_degree, s_pstart, s_edges);
	return degen(s_n, s_peel_sequence, s_core, s_pstart, s_edges, s_degree, s_vis, s_heap,false);
}