#ifndef RSTARTREE_GRAPH_H
#define RSTARTREE_GRAPH_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
typedef unsigned int vid_t;
typedef int num_t;

using namespace std;

class BiGraph
{

public:
    BiGraph(string dir) {
        num_v1 = 0;
        num_v2 = 0;
        num_edges = 0;

        neighbor_v1.clear();
        neighbor_v2.clear();


        degree_v1.clear();
        degree_v2.clear();

        core_v1.clear();
        core_v2.clear();

        //KKCore index left (x,*) right (*,x)
        left_index.clear();
        right_index.clear();
        v1_max_degree = 0;
        v2_max_degree = 0;
        this->dir = dir;
        loadGraph(dir);
    }

    BiGraph();
    ~BiGraph() {}

    void addEdge(vid_t u, vid_t v) {
        neighbor_v1[u].push_back(v);
        ++degree_v1[u];
        if (degree_v1[u] > v1_max_degree) v1_max_degree = degree_v1[u];
        neighbor_v2[v].push_back(u);
        ++degree_v2[v];
        if (degree_v2[v] > v2_max_degree) v2_max_degree = degree_v2[v];
//	num_edges++;
    }
    void deleteEdge(vid_t u, vid_t v);
    bool isEdge(vid_t u, vid_t v);
    num_t getV1Num() { return num_v1; }
    num_t getV2Num() { return num_v2; }
    num_t getV1Degree(vid_t u) { return degree_v1[u]; }
    num_t getV2Degree(vid_t u) { return degree_v2[u]; }
    std::vector<vid_t> & getV2Neighbors(vid_t u) { return neighbor_v2[u]; }
    std::vector<vid_t> & getV1Neighbors(vid_t u) { return neighbor_v1[u]; }
    void print();
    void print(bool hash);
    void printSum();
    void printCout();

public:

    void init(unsigned int num1, unsigned int num2) {
        num_v1 = num1;
        num_v2 = num2;
        num_edges = 0;

        neighbor_v1.resize(num_v1);
        neighbor_v2.resize(num_v2);

        degree_v1.resize(num_v1);
        degree_v2.resize(num_v2);

        fill_n(degree_v1.begin(), num_v1, 0);
        fill_n(degree_v2.begin(), num_v2, 0);

        left_delete.resize(num_v1);
        right_delete.resize(num_v2);
    }

    void loadGraph(string dir) {
        unsigned int n1, n2;
        unsigned int edges = 0;
        int u, v;
        int r;

        string metaFile = dir + ".meta";
        string edgeFile = dir + ".e";

        FILE *metaGraph = fopen(metaFile.c_str(), "r");
        FILE *edgeGraph = fopen(edgeFile.c_str(), "r");

        // load the number of vertices in BiGraph
        if (fscanf(metaGraph, "%d\n%d", &n1, &n2) != 2) {
            fprintf(stderr, "Bad file format: n1 n2 incorrect\n");
            exit(1);
        }

        fprintf(stdout, "n1: %d, n2: %d\n", n1, n2);

        init(n1, n2);

        while ((r = fscanf(edgeGraph, "%d %d", &u, &v)) != EOF) {
            //fprintf(stderr, "%d, %d\n", u, v);
            if (r != 2) {
                fprintf(stderr, "Bad file format: u v incorrect\n");
                exit(1);
            }

            addEdge(u, v);
            num_edges++;
        }

        fclose(metaGraph);
        fclose(edgeGraph);

        for (int i = 0; i < num_v1; ++i) {
            neighbor_v1[i].shrink_to_fit();
            sort(neighbor_v1[i].begin(), neighbor_v1[i].end());

        }
        for (int i = 0; i < num_v2; ++i) {
            neighbor_v2[i].shrink_to_fit();
            sort(neighbor_v2[i].begin(), neighbor_v2[i].end());
        }
    }

    void compressGraph();

    std::string dir;
    num_t num_v1;
    num_t num_v2;
    num_t num_edges;

    std::vector<std::vector<vid_t>> neighbor_v1;
    std::vector<std::vector<vid_t>> neighbor_v2;

    std::vector<std::unordered_set<vid_t>> neighborHash_v1;
    std::vector<std::unordered_set<vid_t>> neighborHash_v2;

    std::vector<int> degree_v1;
    std::vector<int> degree_v2;

    std::vector<num_t> core_v1;
    std::vector<num_t> core_v2;

public:

    //KKCore index left (x,*) right (*,x)
    std::vector<std::vector<int>> left_index;
    std::vector<std::vector<int>> right_index;
    int v1_max_degree;
    int v2_max_degree;
    std::vector<bool> left_delete;
    std::vector<bool> right_delete;
    // for dynamic update
    std::vector<std::vector<int>> left_index_old;
    std::vector<std::vector<int>> right_index_old;
    //BiGraph operator=(const BiGraph& g);

public:
    int get_left_index_with_fixed_left_k(vid_t u, int left_k);
    //BiGraph& operator=(const BiGraph& g_);
};

class subBiGraph {
public:
    std::vector<int> leftVertex;
    std::vector<int> rightVertex;
    subBiGraph(){};
};

#endif //RSTARTREE_GRAPH_H
