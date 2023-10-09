#ifndef RSTARTREE_KLSKYLINECORE_H
#define RSTARTREE_KLSKYLINECORE_H
#include <vector>
#include "DiGraph.h"

class KlSkylineCore {
public:
    int id;
    int k;
    int l;
    std::vector<int> vertices;
    std::vector<vector<int>> subGraphs;
    std::vector<int> addedVertices;
    std::vector<int> edgeVertices;
    vector<pair<int, int>> connectEdges;
    //    std::vector<pair<int, int>> criticalEdges;
    KlSkylineCore(){};
};

void loadKlSkylineCores(ifstream &inFile, vector<KlSkylineCore> &cores, int &kMax, int &lMax) {
    kMax = 0;
    lMax = 0;
    int num;
    inFile >> num;
    cores.resize(num);
    for (int i = 0; i < num; i++) {
        int tempA, tempB;
        inFile >> tempA >> tempB;
        cores[i].id = i;
        cores[i].k = tempA;
        cores[i].l = tempB;
        kMax = max(kMax, tempA);
        lMax = max(lMax, tempB);
        int vertexNum;
        inFile >> vertexNum;
        cores[i].vertices.resize(vertexNum);
        for (int j = 0; j < vertexNum; j++) {
            int vid;
            inFile >> vid;
            cores[i].vertices[j] = vid;
        }
    }
};

#endif //RSTARTREE_KLSKYLINECORE_H
