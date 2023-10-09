#ifndef RSTARTREE_MULTILAYER_GRAPH_H
#define RSTARTREE_MULTILAYER_GRAPH_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <set>


class MultilayerGraph {

public:
    MultilayerGraph();
    MultilayerGraph(std::string dataset_path);
    void load_dataset(std::string dataset_path);
    void add_edge(int from_node, int to_node, int layer);
    std::unordered_map<int, int> order_layers(std::ifstream& dataset_file);
    std::vector<int> get_nodes();
    double get_number_of_edges(int layer = -1);
    std::unordered_map<int, double> get_number_of_edges_layer_by_layer();
    std::vector<std::vector<int>> neighbors;
    int get_layer_mapping(int layer);

public:
    int number_of_layers;
    std::vector<int> layers_iterator;
    std::unordered_map<int, int> layers_map;

    int number_of_nodes;
    int maximum_node;
    std::vector<int> nodes_iterator;
    std::vector<std::vector<std::vector<int>>> adjacency_list;

    std::string dataset_path;
};

class SubGraph {
public:
    std::vector<int> vertices;
    SubGraph(){};
};

class MultiSkylineCore {
public:
    int id;
    std::vector<int> parameters;
    std::vector<int> vertices;
    std::vector<std::vector<int>> subGraphs;
//    std::vector<std::pair<int, int>> criticalEdges;
    std::vector<std::pair<int, int>> connectEdges;
    std::vector<int> addedVertices;
public:
    MultiSkylineCore(){};
};


#endif //RSTARTREE_MULTILAYER_GRAPH_H
