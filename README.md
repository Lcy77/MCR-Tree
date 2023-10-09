# MCR-Tree

## Introduction
This repository contains the code used in our paper: **MCR-Tree: An Efficient Index for Multi-dimensional Core Search**.
MCR-Tree projects graph vertices into a multi-dimensional space by leveraging the skyline property of cores, and  integrates skyline coreness and connectivity information into an R-tree to facilitate the multi-dimensional core search..

## Datasets
We use 12 publicly available real-world networks, their types include directed graphs, bipartite graphs, and multilayer graphs.
Directed graphs are obtained from [SNAP](http://snap.stanford.edu/) and [LAW](https://law.di.unimi.it/datasets.php). 
Bipartite graphs are obtained from [KONECT](http://konect.uni-koblenz.de/networks).
Multilayer graphs are obtained from [Manlio De Domenico's Website](http://deim.urv.cat/~manlio.dedomenico/data.php).

An example format of the input data is shown in folder kl-core/test.

## Algorithms
The following files are the codes for our proposed algorithms.
Due to the inconsistent data format of different types of graphs, we split the code of the MCR-Tree into several different files.
We implemented all codes the using C++ with CLion 2022.1.1.
1. **RStarTree.h** : Algorithms to build R*-Tree [1].
2. **klSearch.cpp** : Algorithms to construct MCR-Tree for the (*k*, *l*)-core and perform (*k*, *l*)-core saerch.
3. **abSearch.cpp** : Algorithms to construct MCR-Tree for the (*&alpha;*, *&beta;*)-core and perform (*&alpha;*, *&beta;*)-core saerch.
4. **multilayerSearch.cpp** : Algorithms to construct MCR-Tree for the multilayer **k**-core and perform multilayer **k**-core saerch.
5. **klUpdate.cpp**, **abUpdate.cpp**, **multilayerUpdate.cpp** : Index maintenance algorithm of MCR-Tree.

[1] Beckmann N, Seeger B. A revised R*-tree in comparison with related index structures[C]//Proceedings of the 2009 ACM SIGMOD International Conference on Management of data. 2009: 799-812.

## Usage
1. **Get skyline corenesses of the graph.** Existing works [2] [3] can be used to perform core decomposition and get skylien corenesses for each vertex. 

In the folder _kl-core/test_, _test_degree.dat_ and _test_graph.dat_ are the original graph, and _test_klMax_ and _test_skyline_ are files related to the skyline coreness.

2. **Build MCR-Tree and perform core search using Algorithms 1~3.** There are two steps to build MCR-Tree: 
(i) build R-tree using skyline corenesses, (ii) compute contents of nodes for MCR-Tree.

For example, after the project is compiled, just type in:
```
/m-core/cmake-build-release/klSearch test
```
MCR-Tree of the test dataset is built and some core searches are performed.

[2] Fang Y, Wang Z, Cheng R, et al. Effective and efficient community search over large directed graphs[J]. IEEE Transactions on Knowledge and Data Engineering, 2018, 31(11): 2093-2107.

[3] Galimberti E, Bonchi F, Gullo F, et al. Core decomposition in multilayer networks: Theory, algorithms, and applications[J]. ACM Transactions on Knowledge Discovery from Data, 2020, 14(1): 1-40.

## Requirements
+ cmake
+ g++
+ OpenMP