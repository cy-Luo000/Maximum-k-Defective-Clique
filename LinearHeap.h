#ifndef _LINEAR_HEAP_H_
#define _LINEAR_HEAP_H_

#include "Utility.h"

class ListLinearHeap {
public:
	int n; // number vertices
	int key_cap; // the maximum allowed key value

	int min_key; // possible min key
	int max_key; // possible max key

	int *key_s; // key of vertices

	int *head_s; // head of doubly-linked list for a specific weight

	int *pre_s; // pre for doubly-linked list
	int *next_s; // next for doubly-linked list

public:
	ListLinearHeap(int _n, int _key_cap) {
		n = _n;
		key_cap = _key_cap;

		min_key = max_key = key_cap;

		head_s = key_s = pre_s = next_s = nullptr;
	}

	~ListLinearHeap() {
		if(head_s != nullptr) {
			delete[] head_s;
			head_s = nullptr;
		}
		if(pre_s != nullptr) {
			delete[] pre_s;
			pre_s = nullptr;
		}
		if(next_s != nullptr) {
			delete[] next_s;
			next_s = nullptr;
		}
		if(key_s != nullptr) {
			delete[] key_s;
			key_s = nullptr;
		}
	}

	void init(int _n, int _key_cap, int *_id_s, int *_key_s) {
		if(key_s == nullptr) key_s = new int[n];
		if(pre_s == nullptr) pre_s = new int[n];
		if(next_s == nullptr) next_s = new int[n];
		if(head_s == nullptr) head_s = new int[key_cap+1];

		//assert(_key_cap <= key_cap);
		min_key = max_key = _key_cap;
		for(int i = 0;i <= _key_cap;i ++) head_s[i] = n;

		for(int i = 0;i < _n;i ++) {
			int id = _id_s[i];
			int key = _key_s[id];
			//assert(id < n); assert(key <= _key_cap);

			key_s[id] = key; pre_s[id] = n; next_s[id] = head_s[key];
			if(head_s[key] != n) pre_s[head_s[key]] = id;
			head_s[key] = id;

			if(key < min_key) min_key = key;
		}
	}

	int get_key(int id) { return key_s[id]; }

	void get_ids(int *vs, int &vs_size) {
		for(int i = min_key;i <= max_key;i ++) {
			for(int id = head_s[i];id != n;id = next_s[id]) {
				vs[vs_size ++] = id;
			}
		}
	}

	bool get_min(int &id, int &key) {// return true if success, return false otherwise
		while(min_key <= max_key&&head_s[min_key] == n) ++ min_key;
		if(min_key > max_key) return false;

		id = head_s[min_key];
		key = min_key;

		//assert(key_s[id] == key);

		return true;
	}

	bool pop_min(int &id, int &key) {// return true if success, return false otherwise
		while(min_key <= max_key&&head_s[min_key] == n) ++ min_key;
		if(min_key > max_key) return false;

		id = head_s[min_key];
		key = min_key;

		key_s[id] = key_cap+1;
		//assert(key_s[id] == key);

		head_s[min_key] = next_s[id];
		if(head_s[min_key] != n) pre_s[head_s[min_key]] = n;
		return true;
	}

	int decrement(int id, int dec) {
		//assert(key_s[id] >= dec);
		if(key_s[id] > key_cap) return 0;

		if(pre_s[id] == n) {
			//assert(head_s[key_s[id]] == id);
			head_s[key_s[id]] = next_s[id];
			if(next_s[id] != n) pre_s[next_s[id]] = n;
		}
		else {
			int pid = pre_s[id];
			next_s[pid] = next_s[id];
			if(next_s[id] != n) pre_s[next_s[id]] = pid;
		}

		int &key = key_s[id];
		key -= dec; pre_s[id] = n; next_s[id] = head_s[key];
		if(head_s[key] != n) pre_s[head_s[key]] = id;
		head_s[key] = id;

		if(key < min_key) min_key = key;
		return key;
	}

	void del(int id) {
		if(key_s[id] > key_cap) return;

		if(pre_s[id] == n) {
			//assert(head_s[key_s[id]] == id);
			head_s[key_s[id]] = next_s[id];
			if(next_s[id] != n) pre_s[next_s[id]] = n;
		}
		else {
			int pid = pre_s[id];
			next_s[pid] = next_s[id];
			if(next_s[id] != n) pre_s[next_s[id]] = pid;
		}
	}
};

#endif