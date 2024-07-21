#ifndef _COLOR_PACKER_H_
#define _COLOR_PACKER_H_

#include "Utility.h"
using namespace std;
class colorPacker
{
private:
	/* data */
public:
	int n;// the length of the color vector
	int *colLab;//the color of the element in color vector(color lable)
	vector<vector<int>> AdjList;//the graph of non mutual pairs
	vector<int> vChose;
	bool *colMtx;
	int colorNum;
	int maxWsum;

	colorPacker(std::vector<std::vector<bool>>& MuEx, int *colVec, int ColSz);
	colorPacker(int n);
	void set(std::vector<std::vector<bool>>& MuEx, int *colVec, int ColSz);
	void reset();
	void clear();
	void coloring();
	void coloring(int *colVec, int* neiInP, int P_end);
	int colPack(int *colVec, int colSz ,int *neiInP,int P_end,int threshold);
	~colorPacker();
};
colorPacker::colorPacker(int n_){// n_ is the upperbound of n
	colLab=new int[n_];
	memset(colLab, 0, n_*sizeof(int));
	AdjList.resize(n_);
	colMtx=new bool[n_*n_];
	memset(colMtx, false, n_*n_*sizeof(bool));
	n=0, colorNum=0, maxWsum=0;
}
void colorPacker::set(std::vector<std::vector<bool>>& MuEx, int *colVec, int ColSz){
	n=ColSz;// n is the color size
	maxWsum=0;
	for (int i = 0; i < ColSz; i++){
		int u=colVec[i];
		for (int j = i+1; j < ColSz; j++){
			int v=colVec[j];
#ifdef _DBUG_
			assert(u<n_ && u>=0);
			assert(v<n_ && v>=0);
#endif
			// printf("u: %d, v: %d, muex[u][v]: %d\n", u,v,MuEx[u][v]);
			if(!MuEx[u][v]) AdjList[i].push_back(j);//if (u,v) is not mutual exclusive, then build the edge
		}
	}
}
void colorPacker::reset(){
	//reset the color lable, adjlist and color mtx
	memset(colLab, 0, n*sizeof(int));
	for (int i = 0; i < n; i++) AdjList[i].clear();
	memset(colMtx, false, n*n*sizeof(bool));
	vChose.clear();
	colorNum=0, n=0, maxWsum=0;
}
colorPacker::colorPacker(std::vector<std::vector<bool>>& MuEx, int *colVec, int ColSz)
{
	for (int i = 0; i < ColSz; i++){
		int u=colVec[i];
		if(u<0){
			printf("u: %d, i: %d, sz: %d\n",u, i, ColSz);
			exit(0);
		}
#ifdef _DBUG_
		assert(u>=0);
		assert(u<n_);
#endif
	}
	
	n=ColSz;// n is the color size
	colorNum=0;
	AdjList.resize(n);
	//construct the graph
	for (int i = 0; i < ColSz; i++){
		int u=colVec[i];
		for (int j = i+1; j < ColSz; j++){
			int v=colVec[j];
#ifdef _DBUG_
			assert(u<n_ && u>=0);
			assert(v<n_ && v>=0);
#endif
			// printf("u: %d, v: %d, muex[u][v]: %d\n", u,v,MuEx[u][v]);
			if(!MuEx[u][v]) AdjList[i].push_back(j);//if (u,v) is not mutual exclusive, then build the edge
		}
	}
	colMtx=new bool[n*n];
	memset(colMtx,false, n*n*sizeof(bool));
	colLab=new int[n];
	memset(colLab,0, n*sizeof(int));
}
void colorPacker::clear(){
	n=0;
	if(colLab!=NULL){
		delete[] colLab;
		colLab=NULL;
	}
	if(AdjList.size()>0){
		AdjList.clear();
	}
	if(!vChose.empty())
		vChose.clear();
	if(colMtx!=NULL){
		delete[] colMtx;
		colMtx=NULL;
	}
}
void colorPacker::coloring(){
	int maxCol=-1;
	// do the coloring

	for (int i = 0; i < n; i++){
		int col=0;
		while (colMtx[n*i+col]) col++;
		if(col>maxCol) {
			maxCol=col;
			vChose.push_back(i);
		}
		colLab[i]=col;
		for (auto j:AdjList[i]) colMtx[n*j+col]=true;
		// vChose.push_back(i);
	}
	colorNum=maxCol+1;
	return;
}
void colorPacker::coloring(int* colVec, int *neiInP, int P_end){
	int maxCol=-1;
	// do the coloring
	for (int i = 0; i < n; i++){
		int col=0;
		while (colMtx[n*i+col]) col++;
		if(col>maxCol) {
			maxCol=col;
			vChose.push_back(i);
			maxWsum+=(P_end-neiInP[colVec[i]]+vChose.size()-1);
		}
		colLab[i]=col;
		for (auto j:AdjList[i]) colMtx[n*j+col]=true;
	}
	colorNum=maxCol+1;
	return;
}
int colorPacker::colPack(int *colVec, int colSz ,int *neiInP,int P_end,int threshold){
	int choseNum=0,mssEdge=0;
	if(threshold>=maxWsum) return vChose.size();
	for (int i = 0; i < vChose.size(); i++){
		int u=colVec[vChose[i]];
		if(mssEdge+i+(P_end-neiInP[u])<=threshold) choseNum++;
		else break;
		mssEdge+=(i+P_end-neiInP[u]);
	}
	return choseNum;
}
colorPacker::~colorPacker()
{
	n=0;
	if(colLab!=NULL){
		delete[] colLab;
		colLab=NULL;
	}
	if(AdjList.size()>0){
		AdjList.clear();
	}
	if(!vChose.empty()) vChose.clear();
	if(colMtx!=NULL){
		delete[] colMtx;
		colMtx=NULL;
	}
}

#endif