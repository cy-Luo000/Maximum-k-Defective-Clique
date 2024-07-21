# Towards A Faster Branching Algorithm for the Maximum k-Defective Clique Problem

Note that: 

In the implementation, there is no need to explicitly build complement graphs , as we use the adjacency matrix to store subgraphs.

The lower bound of DnBk are simply set to at least k+2 to screen out trivial cases.

## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "DnBk", which actually corresponds to the our algorithm for optimization problem (Ours).

## Run the code

```sh
$ ./DnBk {path_to_binary_compressed_graph} {k_value}
```

An example of computing the exact maximum 1-defective clique for the dataset soc-slashdot is as follows
```sh
$ ./DnBk data/soc-slashdot.bin 1
```

The solution is in output file "KDC.txt".

## Data format
We adopt the time-efficient binary format rather than the text format.  Several demonstration graphs are given in "data" folder.

Transforming [network-repo graphs](http://lcs.ios.ac.cn/~caisw/Resource/realworld%20graphs.tar.gz) from text format to binary format is as follows:
```sh
$ g++ ./toBin.cpp -o toBin
$ mv ./data/soc-slashdot ./data/soc-slashdot.clq
$ ./toBin data/soc-slashdot.clq
$ ls data/soc-slashdot.bin
```# Maximum-k-Defective-Clique
