#ifndef RSTARTREE_DIGRAPH_H
#define RSTARTREE_DIGRAPH_H

#include <string>
#include <vector>
#include <fstream>

class DiGraph {
public:
    std::string degree_file;
    std::string graph_file;
    //used in reading data, to create an array with exact length of neighbor's number
    int numOfNodes, numOfEdges;
    std::vector<std::vector<int>> inGraph, outGraph;
    std::vector<int> inNum, outNum;
    std::ifstream graphReader;

    /**
* description: take $sampleRate$ percent vertices from the data
 * when sampling, we need do a reflect to make vertices ID begins from 0 and keep continuous
*/
public:
    DiGraph(const std::string& degreeFile, const std::string& graphFile);
    DiGraph(const std::string &graphFile);
    void read_graph();
    std::vector<std::vector<int>> get_inGraph();
    std::vector<std::vector<int>> get_outGraph();
};

class SubDiGraph {
public:
    std::vector<int> vertices;
};

#endif //RSTARTREE_DIGRAPH_H
