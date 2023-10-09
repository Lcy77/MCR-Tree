#include "DiGraph.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

DiGraph::DiGraph(const std::string& degreeFile, const std::string& graphFile) :
        degree_file(degreeFile), graph_file(graphFile), numOfNodes(0), numOfEdges(0) {
    try {
        graphReader.open(graph_file);
        std::string tempString;
        std::getline(graphReader, tempString);
        std::stringstream ss(tempString);
        std::string token;
        std::getline(ss, token, ' ');
        this->numOfNodes = std::stoi(token);
        inGraph.resize(numOfNodes);
        outGraph.resize(numOfNodes);
        inNum.resize(numOfNodes);
        outNum.resize(numOfNodes);
    } catch (const std::ifstream::failure& e) {
        std::cerr << "Error: Failed to open file." << std::endl;
        return;
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: Invalid argument." << std::endl;
        return;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: Out of range." << std::endl;
        return;
    } catch (...) {
        std::cerr << "Error: Unknown exception." << std::endl;
        return;
    }
}


//fast read file to load the outGraph matrix
DiGraph::DiGraph(const std::string &graphFile) {
    std::ifstream graphReader(graphFile);
    if (graphReader.is_open()) {
        std::string tempString;
        std::getline(graphReader, tempString);
        std::istringstream iss(tempString);
        iss >> numOfNodes;

        outGraph.resize(numOfNodes);
        outNum.resize(numOfNodes);

        std::string tempString2;
        while (std::getline(graphReader, tempString2)) {
            std::istringstream iss2(tempString2);
            int nid;
            iss2 >> nid;
            outGraph[nid].resize(tempString2.length() - 1);

            int neighborId;
            int k = 0;
            while (iss2 >> neighborId) {
                if (neighborId == nid) {
                    continue;
                }
                outGraph[nid][outNum[nid]++] = neighborId;
                ++k;
            }
        }

        graphReader.close();
    }
    else {
        std::cerr << "Failed to open graph file: " << graphFile << std::endl;
    }
}

std::vector<std::string> split_string(const std::string& input) {
    std::vector<std::string> tokens;
    std::istringstream iss(input);
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

void DiGraph::read_graph() {
    std::ifstream edgeReader(degree_file);
    if (edgeReader.is_open()) {
        std::string tempString1, tempString2;
        std::vector<std::string> split1, split2;
        while (std::getline(edgeReader, tempString1)) {
            split1 = split_string(tempString1);
            int nid = std::stoi(split1[0]);
            inGraph[nid].resize(std::stoi(split1[1]));
            outGraph[nid].resize(std::stoi(split1[2]));
        }
        edgeReader.close();

        if (graphReader.is_open()) {
            while (std::getline(graphReader, tempString2)) {
                split2 = split_string(tempString2);
                int nid = std::stoi(split2[0]);
                for (int k = 1; k < split2.size(); k++) {
                    int neighborId = std::stoi(split2[k]);
                    if (neighborId == nid) {
                        continue;
                    }
                    outGraph[nid][outNum[nid]++] = neighborId;
                    inGraph[neighborId][inNum[neighborId]++] = nid;
                }
                outGraph[nid].resize(outNum[nid]);
            }
            graphReader.close();
        }
    }
    else {
        std::cerr << "Failed to open degree file: " << degree_file << std::endl;
    }
}


std::vector<std::vector<int>> DiGraph::get_inGraph() {
    return inGraph;
}

std::vector<std::vector<int>> DiGraph::get_outGraph() {
    return outGraph;
}