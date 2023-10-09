#include "multilayer_graph.h"

#include "multilayer_graph.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

MultilayerGraph::MultilayerGraph() {
    number_of_layers = 0;
    layers_iterator.resize(0);
    layers_map = {};

    number_of_nodes = 0;
    maximum_node = 0;
    nodes_iterator.resize(0);
    adjacency_list = {};
}

MultilayerGraph::MultilayerGraph(std::string dataset_path) {
    number_of_layers = 0;
    layers_iterator.resize(0);
    layers_map = {};

    number_of_nodes = 0;
    maximum_node = 0;
    nodes_iterator.resize(0);
    adjacency_list = {};

    load_dataset(dataset_path);
    this->dataset_path = dataset_path;
}

void MultilayerGraph::load_dataset(std::string dataset_path) {
    std::ifstream dataset_file;

    dataset_file.open(dataset_path);

    std::string first_line;
    std::getline(dataset_file, first_line);
    std::vector<std::string> split_first_line;
    std::istringstream ss(first_line);
    std::string token;
    while (std::getline(ss, token, ' ')) {
        split_first_line.push_back(token);
    }

    number_of_layers = std::stoi(split_first_line[0]);
    layers_iterator.resize(number_of_layers);
    std::iota(layers_iterator.begin(), layers_iterator.end(), 0);

    number_of_nodes = std::stoi(split_first_line[1]);
    maximum_node = std::stoi(split_first_line[2]);
    nodes_iterator.resize(maximum_node + 1);
    std::iota(nodes_iterator.begin(), nodes_iterator.end(), 0);

    adjacency_list.resize(maximum_node + 1);
    for (auto &neighbors: adjacency_list) {
        neighbors.resize(number_of_layers);
    }

    auto layers_map = order_layers(dataset_file);

    std::string line;
    while (std::getline(dataset_file, line)) {
        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> split_line;
        while (std::getline(ss, token, ' ')) {
            split_line.push_back(token);
        }

        int layer = std::stoi(split_line[0]);
        int from_node = std::stoi(split_line[1]);
        int to_node = std::stoi(split_line[2]);

        add_edge(from_node, to_node, layers_map[layer]);
    }
    dataset_file.close();
}

void MultilayerGraph::add_edge(int from_node, int to_node, int layer) {
    if (from_node != to_node) {
        adjacency_list[from_node][layer].push_back(to_node);
        adjacency_list[to_node][layer].push_back(from_node);
    }
}

std::unordered_map<int, int> MultilayerGraph::order_layers(std::ifstream &dataset_file) {
    // map of the layers
    std::unordered_map<int, int> layers_map;
    int layers_oracle = 0;
    std::string line;
    while (std::getline(dataset_file, line)) {
        std::istringstream ss(line);
        std::string token;
        std::getline(ss, token, ' ');
        int layer = std::stoi(token);

        if (layers_map.find(layer) == layers_map.end()) {
            layers_map[layer] = layers_oracle;
            this->layers_map[layers_oracle] = layer;
            layers_oracle++;
        }
    }
    dataset_file.clear();
    dataset_file.seekg(0, std::ios::beg);
    std::getline(dataset_file, line);
    return layers_map;
}

std::vector<int> MultilayerGraph::get_nodes() {
    std::vector<int> nodes;
    if (number_of_nodes == maximum_node) {
        for (int i = 1; i <= maximum_node; i++) {
            nodes.push_back(i);
        }
    } else {
        for (int i = 0; i <= maximum_node; i++) {
            nodes.push_back(i);
        }
    }
    return nodes;
}

double MultilayerGraph::get_number_of_edges(int layer) {
    double number_of_edges = 0.0;
    for (const auto &neighbors: adjacency_list) {
        for (int inner_layer = 0; inner_layer < number_of_layers; inner_layer++) {
            if (layer == -1 || layer == inner_layer) {
                number_of_edges += neighbors[inner_layer].size();
            }
        }
    }

    return number_of_edges / 2;
}

std::unordered_map<int, double> MultilayerGraph::get_number_of_edges_layer_by_layer() {
    std::unordered_map<int, double> number_of_edges_layer_by_layer;
    for (int layer: layers_iterator) {
        double layer_edge_count = 0.0;
        for (const auto &neighbors: adjacency_list) {
            layer_edge_count += neighbors[layer].size();
        }
        number_of_edges_layer_by_layer[layer] = layer_edge_count / 2;
    }
    return number_of_edges_layer_by_layer;
}

int MultilayerGraph::get_layer_mapping(int layer) {
    return layers_map[layer];
}
