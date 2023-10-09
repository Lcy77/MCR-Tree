#ifndef RSTARTREE_SKYLINECORES_H
#define RSTARTREE_SKYLINECORES_H

#include <vector>
#include <fstream>
#include <unordered_set>
#include "graph.h"
using namespace std;
class biSkylineCore {
public:
    int id;
    int a, b;
    vector<int> vertexU;
    vector<int> vertexV;
    vector<int> vertices;
    vector<int> coreNeighbor;
    vector<int> neighborU;
    vector<int> neighborV;
    vector<int> addedVertices;
    vector<vector<int>> subGraphs;
    vector<pair<int, int>> connectEdges;
    biSkylineCore(){}
};

void loadBiSkylineCores(ifstream &inFile, vector<biSkylineCore> &cores, int &alphaMax, int &betaMax) {
    alphaMax = 0;
    betaMax = 0;
    int num;
    inFile >> num;
    cores.resize(num);
    for (int i = 0; i < num; i++) {
//        printf("%d / %d\n", i, num);
        int tempA, tempB;
        inFile >> tempA >> tempB;
        cores[i].id = i;
        cores[i].a = tempA;
        cores[i].b = tempB;
        alphaMax = max(alphaMax, tempA);
        betaMax = max(betaMax, tempB);
        int numU;
        inFile >> numU;
        cores[i].vertexU.resize(numU);
        for (int j = 0; j < numU; j++) {
            int uid;
            inFile >> uid;
            cores[i].vertexU[j] = uid;
        }
        int numV;
        inFile >> numV;
        cores[i].vertexV.resize(numV);
        for (int j = 0; j < numV; j++) {
            int vid;
            inFile >> vid;
            cores[i].vertexV[j] = vid;
        }
    }
}

#endif //RSTARTREE_SKYLINECORES_H
