#include <string>
#include <ctime>
#include <sys/time.h>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <stdio.h>
#include <queue>
#include <stack>
#include <omp.h>
#include "RStarTree.h"
#include "skylineCores.h"
#include "graph.h"
#include "utilities/UnionFind.h"

using namespace std;
const int dimensions = 2;
const int min_child_items = 32;
const int max_child_items = 64;
typedef RStarTree<int, dimensions, min_child_items, max_child_items> RTree;
typedef RTree::BoundingBox BoundingBox;

int alphaMax, betaMax;
int alpha, beta;
int uNum, vNum;
vector<biSkylineCore> skylineCores;
map<int, vector<int>> vertexCore;
map<pair<int, int>, int> mapCoreId;
vector<int> alphaBetaMax;
map<pair<int, int>, int> mapCoCore;
map<pair<int, int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>> mapCoNode;
vector<vector<int>> coreMatrix;
vector<vector<int>> nodeMatrix;
string dataset;

int pointNum;

BoundingBox bounds(int x, int y, int w, int h) {
    BoundingBox bb;

    bb.edges[0].first = x;
    bb.edges[0].second = x + w;

    bb.edges[1].first = y;
    bb.edges[1].second = y + h;

    return bb;
}

struct UpdatedVertex {
    int vertexID;
    vector<pair<int, int>> newCores;
    vector<pair<int, int>> oldCores;
};

struct addOneEdge {
    int fromVertex;
    int toVertex;
    vector<UpdatedVertex> updatedVertices;
    vector<pair<int, int>> vertexPairs;
};

struct deleteOneEdge {
    int fromVertex;
    int toVertex;
    vector<UpdatedVertex> updatedVertices;
    vector<pair<int, int>> vertexPairs;
};

struct Visitor {
    int count;
    bool ContinueVisiting;

    Visitor() : count(0), ContinueVisiting(true) {};

    void operator()(const RTree::Leaf *const leaf) {
        count++;
    }
};

bool isEqual(std::vector<int> const &a, std::vector<int> const &b) {
    if (a.size() != b.size()) return false;
    return std::equal(a.begin(), a.end(), b.begin());
}

void preProcess(vector<biSkylineCore> &skylineCores, BiGraph &g) {
    for (int i = 0; i < skylineCores.size(); i++) {
        if (i % 10000 == 0) {
            printf("%d / %d\n", i, skylineCores.size());
        }
        mapCoreId.insert(make_pair(make_pair(skylineCores[i].a, skylineCores[i].b), i));
        for (auto u : skylineCores[i].vertexU) {
            vertexCore[u].push_back(i);
        }
        for (auto v : skylineCores[i].vertexV) {
            vertexCore[v + uNum].push_back(i);
        }
    }
}

void loadAlphaBetaMax(ifstream &ifstream1) {
    alphaBetaMax.resize(alphaMax + 1);
    coreMatrix.resize(alphaMax + 1);
    nodeMatrix.resize(alphaMax + 1);
    //// load alpha beta max file
    for (int i = 1; i <= alphaMax; i++) {
        int alpha;
        int beta;
        ifstream1 >> alpha >> beta;
        alphaBetaMax[alpha] = beta - 1;
        pointNum += (beta - 1);
        coreMatrix[alpha].resize(beta);
        nodeMatrix[alpha].resize(beta);
        std::fill(coreMatrix[alpha].begin(), coreMatrix[alpha].end(), -1);
        std::fill(nodeMatrix[alpha].begin(), nodeMatrix[alpha].end(), -1);
    }
}

void rtreeDFS(RStarTree<int, dimensions, min_child_items, max_child_items>::Node * node, int &coreCnt, vector<biSkylineCore> &newSkylineCores) {
    mapCoNode[{node->bound.edges[0].first, node->bound.edges[1].first}].push_back(node);
    nodeMatrix[node->bound.edges[0].first][node->bound.edges[1].first] = 1;
    if (node->hasLeaves) {
        node->left = coreCnt;
        for (auto item1 : node->items) {
            static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Leaf *>(item1)->coreId = coreCnt;
            newSkylineCores.push_back(skylineCores[mapCoreId[{item1->bound.edges[0].first, item1->bound.edges[1].first}]]);
            newSkylineCores.back().id = coreCnt;
            mapCoCore[{item1->bound.edges[0].first, item1->bound.edges[1].first}] = coreCnt;
            coreCnt++;
        }
        node->right = coreCnt-1;
    } else {
        for (auto item1 : node->items) {
            rtreeDFS(static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(item1), coreCnt, newSkylineCores);
        }
        node->left = static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(node->items[0])->left;
        node->right = static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(node->items.back())->right;
//        cout << node->left << " " << node->right << endl;
    }
}

void resortCores(RTree &tree) {
    auto node = tree.m_root;
    int coreCnt = 0;
    vector<biSkylineCore> newSkylineCores;
//    map<pair<int, int>, int> newMapCoreId;
    if (node != NULL) {
        rtreeDFS(node, coreCnt, newSkylineCores);
    }
    skylineCores = newSkylineCores;
    for (int i = 0; i < skylineCores.size(); i++) {
        coreMatrix[skylineCores[i].a][skylineCores[i].b] = i;
    }
}

void updatedRtreeDFS(RStarTree<int, dimensions, min_child_items, max_child_items>::Node *node, int &coreCnt,
                     map<pair<int, int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>> &updatedMapCoNode,
                     vector<vector<int>> &updatedNodeMatrix, map<pair<int, int>, int> &updatedMapCoCore,
                     map<pair<int, int>, int> &newMapCoCore, vector<biSkylineCore> &updatedSkylineCores,
                     vector<biSkylineCore> &newSkylineCores) {
    updatedMapCoNode[{node->bound.edges[0].first, node->bound.edges[1].first}].push_back(node);
    updatedNodeMatrix[node->bound.edges[0].first][node->bound.edges[1].first] = 1;
    if (node->hasLeaves) {
        node->left = coreCnt;
        for (auto item1: node->items) {
            static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Leaf *>(item1)->coreId = coreCnt;
            newSkylineCores.push_back(
                    updatedSkylineCores[updatedMapCoCore[{item1->bound.edges[0].first, item1->bound.edges[1].first}]]);
            newSkylineCores.back().id = coreCnt;
            newMapCoCore[{item1->bound.edges[0].first, item1->bound.edges[1].first}] = coreCnt;
            coreCnt++;
        }
        node->right = coreCnt - 1;
    } else {
        for (auto item1: node->items) {
            updatedRtreeDFS(static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(item1),
                            coreCnt, updatedMapCoNode, updatedNodeMatrix, updatedMapCoCore, newMapCoCore,
                            updatedSkylineCores, newSkylineCores);
        }
        node->left = static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(node->items[0])->left;
        node->right = static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(node->items.back())->right;
//        cout << node->left << " " << node->right << endl;
    }
}

void updatedResortCores(RTree &tree,
                        map<pair<int, int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>> &updatedMapCoNode,
                        vector<vector<int>> &updatedNodeMatrix, map<pair<int, int>, int> &updatedMapCoCore,
                        vector<vector<int>> &updatedCoreMatrix, vector<biSkylineCore> &updatedSkylineCores,
                        map<int, vector<pair<int, int>>> &updatedMapCoreNode) {
    auto node = tree.m_root;
    int coreCnt = 0;
//    updatedMapCoNode.clear();
    map<pair<int, int>, int> newMapCoCore;
    for (int i = 0; i < updatedNodeMatrix.size(); i++) {
        fill_n(updatedNodeMatrix[i].begin(), updatedNodeMatrix[i].size(), -1);
    }
    vector<biSkylineCore> newSkylineCores;
    if (node != NULL) {
        updatedRtreeDFS(node, coreCnt, updatedMapCoNode, updatedNodeMatrix, updatedMapCoCore, newMapCoCore,
                        updatedSkylineCores, newSkylineCores);
    }
    updatedSkylineCores = newSkylineCores;
    updatedMapCoCore = newMapCoCore;
    for (int i = 0; i < updatedCoreMatrix.size(); i++) {
        fill_n(updatedCoreMatrix[i].begin(), updatedNodeMatrix[i].size(), -1);
    }
    for (int i = 0; i < updatedSkylineCores.size(); i++) {
        updatedCoreMatrix[updatedSkylineCores[i].a][updatedSkylineCores[i].b] = i;
    }
    for (auto item: updatedMapCoNode) {
        for (auto treeNode: item.second) {
            for (int i = treeNode->left; i <= treeNode->right; i++) {
                updatedMapCoreNode[i].push_back({treeNode->bound.edges[0].first, treeNode->bound.edges[1].first});
            }
        }
    }
}

bool compareSkylineCore(const biSkylineCore &a, const biSkylineCore &b) {
    if (a.b != b.b) {
        return a.b > b.b;
    } else {
        return a.a < b.a;
    }
}

//bool compareRNodes(RStarTree<int, dimensions, min_child_items, max_child_items>::Node * node1, RStarTree<int, dimensions, min_child_items, max_child_items>::Node * node2) {
//    if (node1->bound.edges[0].second != node2->bound[0].)
//}

void preProcess2(BiGraph &g, RTree &tree) {
    alphaBetaMax.reserve(alphaMax + 1);
    //// get beta max for each alpha
    for (int i = 1; i <= alphaMax; i++) {
        cout << i << "/" << alphaMax << endl;
        BoundingBox bound;

        for (int j = betaMax; j >= 1; j--) {
            vector<int> cores;  // dominate core for the current core
            vector<int> tempResults;
            bound = bounds(i, j, alphaMax, betaMax);
            tree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
            if (cores.size() != 0) {
                alphaBetaMax[i] = j;
                pointNum += j;
                break;
            }
        }
    }
}

// add Tags and Connected Edges
int outEdgeNum = 0;

double RtreeTime = 0;
void connectivity4(BiGraph &g, RTree &tree) {
    timeval start_time, end_time;
    int cnt = 0;
    UnionFind UF(uNum + vNum);
    UnionFind UF2(uNum + vNum);
    vector<bool> nodes(uNum + vNum);
    std::fill_n(nodes.begin(), nodes.size(), false);
    vector<int> preResult;
    for (auto &skylineCore :skylineCores) {
        for (int i = 0; i < skylineCore.vertexU.size(); i++) {
            skylineCore.vertices.push_back(skylineCore.vertexU[i]);
        }
        for (int i = 0; i < skylineCore.vertexV.size(); i++) {
            skylineCore.vertices.push_back(skylineCore.vertexV[i] + uNum);
        }
    }
    for (int a = alphaMax; a > 0; a--) {
        cout << "Alpha:" << a << "/" << alphaMax << endl;
        if (a < alphaMax && alphaBetaMax[a] == alphaBetaMax[a+1]) {
            bool isBlank = true;
            for (int i = 1; i <= alphaBetaMax[a]; i++) {
//                cout << nodeMatrix[a][i] << endl;
                if (coreMatrix[a][i] != -1 || nodeMatrix[a][i] != -1) {
                    isBlank = false;
                    break;
                }
            }
//            cout << isBlank << endl;
            if (isBlank) {
                cnt += alphaBetaMax[a];
                continue;
            }
        }
        UF.Init(preResult);
        UF2.Init(preResult);
        for (auto &element : preResult) {
            nodes[element] = false;
        }
        preResult.clear();
//        BoundingBox bound;
//        vector<pair<int, int>> preCriticalEdges;
        for (int b = alphaBetaMax[a]; b > 0; b--) {
            if (cnt % 10000 == 0)
                cout << cnt << "/" << pointNum << endl;
            cnt++;
            vector<int> cores;  // dominate core for the current core
            vector<int> tempResults;
            gettimeofday(&start_time, nullptr); // record start time
//            bound = bounds(a, b, alphaMax, 0);
//            tree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
            //// hash table
            for (int tempA = a; tempA <= alphaMax; tempA++) {
                if (b > alphaBetaMax[tempA]) {
                    break;
                }
                int coreId = coreMatrix[tempA][b];
                if (coreId != -1) {
                    cores.push_back(coreId);
                }
            }
            gettimeofday(&end_time, nullptr); // record end time
            RtreeTime += end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1e6; // 计算时间间隔
            unordered_map<int, int> newVertexCore;
            //// get the newly added vertices
            for (int j = 0; j < cores.size(); j++) {
                int coreId = cores[j];
                vector<int> &vertex_u = skylineCores[coreId].vertexU;
                vector<int> &vertex_v = skylineCores[coreId].vertexV;
                if (vertex_u.size() > 0) {
                    auto it = vertex_u.begin();
                    while (it != vertex_u.end()) {
                        if (!nodes[*it]) {
                            nodes[*it] = true;
                            preResult.push_back(*it);
                            tempResults.push_back(*it);
                            newVertexCore[*it] = coreId;
                        }
                        ++it;
                    }
                }
                if (vertex_v.size() > 0) {
                    auto it = vertex_v.begin();
                    while (it != vertex_v.end()) {
                        if (!nodes[*it + uNum]) {
                            nodes[*it + uNum] = true;
                            preResult.push_back(*it + uNum);
                            tempResults.push_back(*it + uNum);
                            newVertexCore[*it + uNum] = coreId;
                        }
                        ++it;
                    }
                }
            }

            //// Calculate connectivity
            for (int j = 0; j < tempResults.size(); j++) {
                if (tempResults[j] < uNum) {
                    for (int node: g.neighbor_v1[tempResults[j]]) {
                        node += uNum;
                        if (nodes[node] && UF.Find(tempResults[j]) != UF.Find(node)) {
                            UF.Unite(tempResults[j], node);
                        }
                    }
                } else {
                    for (int node: g.neighbor_v2[tempResults[j] - uNum]) {
                        if (nodes[node] && UF.Find(tempResults[j]) != UF.Find(node)) {
                            UF.Unite(tempResults[j], node);
                        }
                    }
                }
            }

            //// Count connected subgraphs of current core
            if (coreMatrix[a][b] != -1) {
                int coreID = coreMatrix[a][b];
                map<int, vector<int>> rootVec;
                for (auto vertex: skylineCores[coreID].vertexU) {
                    int tempRoot = UF.Find(vertex);
                    rootVec[tempRoot].push_back(vertex);
                }
                for (auto vertex: skylineCores[coreID].vertexV) {
                    int tempRoot = UF.Find(vertex + uNum);
                    rootVec[tempRoot].push_back(vertex + uNum);
                }

                for (auto &item1: rootVec) {
                    skylineCores[coreID].subGraphs.push_back(item1.second);
                }
            }

            //// Count the connected edges
            for (int j = 0; j < cores.size(); j++) {
                int coreId = cores[j];
                for (auto &subGraph : skylineCores[coreId].subGraphs) {
                    auto it = subGraph.begin();
                    int firstVertex = *it;
                    it++;
                    while (it != subGraph.end()) {
                        UF2.Unite(firstVertex, *it);
                        it++;
                    }
                }
                for (auto &edge : skylineCores[coreId].connectEdges) {
                    UF2.Unite(edge.first, edge.second);
                }
            }
            for (int j = 0; j < tempResults.size(); j++) {
                if (tempResults[j] < uNum) {
                    for (int node: g.neighbor_v1[tempResults[j]]) {
                        node += uNum;
                        if (nodes[node] && UF2.Find(tempResults[j]) != UF2.Find(node)) {
                            UF2.Unite(tempResults[j], node);
                            int tempCoreId = newVertexCore[node];
                            skylineCores[tempCoreId].connectEdges.push_back({tempResults[j], node});
                        }
                    }
                } else {
                    for (int node: g.neighbor_v2[tempResults[j] - uNum]) {
                        if (nodes[node] && UF2.Find(tempResults[j]) != UF2.Find(node)) {
                            UF2.Unite(tempResults[j], node);
                            int tempCoreId = newVertexCore[node];
                            skylineCores[tempCoreId].connectEdges.push_back({tempResults[j], node});
                        }
                    }
                }
            }

            //// calculate tags and out edges for R-tree node
            if (mapCoNode.find({a, b}) != mapCoNode.end()) {
                auto RNodes = mapCoNode.find({a, b})->second;
                vector<bool> tempNodes;
                tempNodes.resize(uNum + vNum);
                vector<int> subgraphRoot;
                vector<int> tempTag;
                for (auto RNode : RNodes) {
                    int coreLeftId = RNode->left;
                    int coreRightId = RNode->right;
                    // get vertices in the range of R-Node and get the tag
                    for (int j = coreLeftId; j <= coreRightId; j++) {
                        for (auto &tempSubgraph : skylineCores[j].subGraphs) {
                            for (auto &tempVertex : tempSubgraph) {
                                tempNodes[tempVertex] = true;
                            }
                            int tag = -1;
                            int tempRoot = UF.Find(tempSubgraph[0]);
                            for (int k = 0; k < subgraphRoot.size(); k++) {
                                if (subgraphRoot[k] == tempRoot) {
                                    tag = k;
                                    break;
                                }
                            }
                            if (tag == -1) {
                                subgraphRoot.push_back(tempRoot);
                                tag = subgraphRoot.size() - 1;
                            }
                            tempTag.push_back(tag);
                        }
                    }
                    RNode->tag = tempTag;
                    // get out edges of R-Node
                    for (int j = coreLeftId; j <= coreRightId; j++) {
                        for (auto &tempEdge : skylineCores[j].connectEdges) {
                            if (!tempNodes[tempEdge.first] || !tempNodes[tempEdge.second]) {
                                RNode->outEdges.push_back(tempEdge);
                                outEdgeNum++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void addOne(RTree tree, BiGraph &g, int addNum) {
    //// load the file
    string addOnePath = "./ab-core/" + dataset + "/" + dataset + "-add";
    addOnePath += to_string(addNum);
    ifstream addOneFile(addOnePath);
    string line;
    vector<addOneEdge> testCase;
    while (std::getline(addOneFile, line)) {
        if (line.size() <= 0)
            break;
        if (line[0] == 'A') {
            int pos = line.find(":");
            line = line.substr(pos + 2);
            istringstream iss(line);
            int fromVertex, toVertex;
            addOneEdge tempAdd;
            while (iss >> fromVertex >> toVertex) {
//                cout << fromVertex << " " << toVertex << endl;
                tempAdd.vertexPairs.push_back({fromVertex, toVertex});
            }
            testCase.push_back(tempAdd);
        } else {
            string delimiter1 = ":";
            string delimiter2 = "|";
            int pos1 = line.find(delimiter1);
            int pos2 = line.find(delimiter2);
            string strID = line.substr(0, pos1);
            string strCoordinatesNew = line.substr(pos1 + 1, pos2 - pos1 - 1);
            string strCoordinatesOld = line.substr(pos2 + 1);

            UpdatedVertex tempVertex;
            tempVertex.vertexID = stoi(strID);

            string delimiter = ")(";

            //// Process the new coordinates
            strCoordinatesNew = strCoordinatesNew.substr(1, strCoordinatesNew.length() - 2);
            size_t pos = 0;
            std::string token;
            while ((pos = strCoordinatesNew.find(delimiter)) != std::string::npos) {
                token = strCoordinatesNew.substr(0, pos); // Extract a coordinate pair string
                // Extract the x and y strings
                size_t commaPos = token.find(",");
                std::string x_str = token.substr(0, commaPos);
                std::string y_str = token.substr(commaPos + 1);
                // Convert strings to integers
                int x = std::stoi(x_str);
                int y = std::stoi(y_str);
                tempVertex.newCores.push_back({x, y});
                strCoordinatesNew.erase(0, pos + delimiter.length()); // Remove the processed part
            }
            int commaPos = strCoordinatesNew.find(",");
            std::string x_str = strCoordinatesNew.substr(0, commaPos);
            std::string y_str = strCoordinatesNew.substr(commaPos + 1);
            int x = std::stoi(x_str);
            int y = std::stoi(y_str);
            tempVertex.newCores.push_back({x, y});

            //// Process the old coordinates
            strCoordinatesOld = strCoordinatesOld.substr(1, strCoordinatesOld.length() - 2);
            pos = 0;
            while ((pos = strCoordinatesOld.find(delimiter)) != std::string::npos) {
                token = strCoordinatesOld.substr(0, pos); // Extract a coordinate pair string
                // Extract the x and y strings
                commaPos = token.find(",");
                x_str = token.substr(0, commaPos);
                y_str = token.substr(commaPos + 1);
                // Convert strings to integers
                x = std::stoi(x_str);
                y = std::stoi(y_str);
                tempVertex.oldCores.push_back({x, y});
                strCoordinatesOld.erase(0, pos + delimiter.length()); // Remove the processed part
            }
            commaPos = strCoordinatesOld.find(",");
            x_str = strCoordinatesOld.substr(0, commaPos);
            y_str = strCoordinatesOld.substr(commaPos + 1);
            x = std::stoi(x_str);
            y = std::stoi(y_str);
            tempVertex.oldCores.push_back({x, y});
            testCase[testCase.size() - 1].updatedVertices.push_back(tempVertex);
        }
    }
    addOneFile.close();

    int cnt = 1;
    clock_t start_time, end_time;
    double total_time = 0;
    for (auto tempCase: testCase) {
        cout << cnt << "/" << testCase.size() << endl;
        cnt++;
        // Get skyline cores in updated vertices
        auto &updatedVertices = tempCase.updatedVertices;
        if (updatedVertices.size() == 0)
            continue;
        //// build temp index for test
        auto updatedSkylineCores = skylineCores;
        map<pair<int, int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>> updatedMapCoNode;
        map<pair<int, int>, int> updatedMapCoCore = mapCoCore;
        vector<vector<int>> updatedCoreMatrix = coreMatrix;
        vector<vector<int>> updatedNodeMatrix = nodeMatrix;
        RTree originTree;
        for (int i = 0; i < updatedSkylineCores.size(); i++) {
            int xx = updatedSkylineCores[i].a;
            int yy = updatedSkylineCores[i].b;
            originTree.Insert(i, bounds(xx, yy, 0, 0));
        }

        start_time = clock();
        map<pair<int, int>, vector<int>> mapOldSkylineCore;
        map<pair<int, int>, vector<int>> mapNewSkylineCore;
        for (auto &tempVertex: updatedVertices) {
            for (int i = 0; i < tempVertex.oldCores.size(); i++) {
                mapOldSkylineCore[{tempVertex.oldCores[i].first, tempVertex.oldCores[i].second}].push_back(
                        tempVertex.vertexID);
            }
            for (int i = 0; i < tempVertex.newCores.size(); i++) {
                mapNewSkylineCore[{tempVertex.newCores[i].first, tempVertex.newCores[i].second}].push_back(
                        tempVertex.vertexID);
            }
        }
        // build 2 R-tree for updated vertices
        vector<biSkylineCore> oldSkylineCore;
        vector<biSkylineCore> newSkylineCore;
        RTree oldTree;
        RTree newTree;
        map<pair<int, int>, int> mapOldCoreId;
        map<pair<int, int>, int> mapNewCoreId;
        int index1 = 0;
        int index2 = 0;
        for (auto &core: mapOldSkylineCore) {
            int xx = core.first.first;  // k
            int yy = core.first.second; // l
            oldTree.Insert(index1, bounds(xx, yy, 0, 0));
            biSkylineCore tempSkylineCore = biSkylineCore();
            tempSkylineCore.a = core.first.first;
            tempSkylineCore.b = core.first.second;
            tempSkylineCore.vertices = core.second;
            oldSkylineCore.push_back(tempSkylineCore);
            mapOldCoreId[{tempSkylineCore.a, tempSkylineCore.b}] = index1;
            index1++;
        }
        for (auto &core: mapNewSkylineCore) {
            int xx = core.first.first;  // k
            int yy = core.first.second; // l
            newTree.Insert(index2, bounds(xx, yy, 0, 0));
            biSkylineCore tempSkylineCore = biSkylineCore();
            tempSkylineCore.a = core.first.first;
            tempSkylineCore.b = core.first.second;
            tempSkylineCore.vertices = core.second;
            newSkylineCore.push_back(tempSkylineCore);
            mapNewCoreId[{tempSkylineCore.a, tempSkylineCore.b}] = index2;
            index2++;
        }
        long long size = 0;
        oldTree.indexSize(RTree::AcceptAny(), Visitor(), size, mapOldCoreId);
        newTree.indexSize(RTree::AcceptAny(), Visitor(), size, mapNewCoreId);

        // updated original skyline core
        for (auto vertex: updatedVertices) {
            // find the changed skyline cores of vertex
            std::sort(vertex.newCores.begin(), vertex.newCores.end());
            std::sort(vertex.oldCores.begin(), vertex.oldCores.end());
            vector<pair<int, int>> coresAdded;
            vector<pair<int, int>> coresDeleted;
            std::set_difference(vertex.newCores.begin(), vertex.newCores.end(),
                                vertex.oldCores.begin(), vertex.oldCores.end(),
                                std::back_inserter(coresAdded));
            std::set_difference(vertex.oldCores.begin(), vertex.oldCores.end(),
                                vertex.newCores.begin(), vertex.newCores.end(),
                                std::back_inserter(coresDeleted));
            // Delete vertices in skyline cores that no longer exist
            for (int i = 0; i < coresDeleted.size(); i++) {
                int coreId = updatedMapCoCore[coresDeleted[i]];
                // remove the vertex in the old core
                auto oldEnd = std::remove(updatedSkylineCores[coreId].vertices.begin(),
                                          updatedSkylineCores[coreId].vertices.end(),
                                          vertex.vertexID);
                updatedSkylineCores[coreId].vertices.erase(oldEnd, updatedSkylineCores[coreId].vertices.end());
                if (updatedSkylineCores[coreId].vertices.size() == 0) {
                    BoundingBox bound = bounds(coresDeleted[i].first, coresDeleted[i].second, 0, 0);
                    originTree.RemoveBoundedArea(bound);
                    updatedMapCoCore.erase({coresDeleted[i].first, coresDeleted[i].second});
                    updatedCoreMatrix[coresDeleted[i].first][coresDeleted[i].second] = -1;
                }
                for (int j = 0; j < updatedSkylineCores[coreId].subGraphs.size(); j++) {
                    auto iter = std::find(updatedSkylineCores[coreId].subGraphs[j].begin(),
                                          updatedSkylineCores[coreId].subGraphs[j].end(), vertex.vertexID);
                    if (iter != updatedSkylineCores[coreId].subGraphs[j].end()) {
                        updatedSkylineCores[coreId].subGraphs[j].erase(iter);
                        if (updatedSkylineCores[coreId].subGraphs[j].size() == 0) {
                            updatedSkylineCores[coreId].subGraphs.erase(
                                    updatedSkylineCores[coreId].subGraphs.begin() + j);
                        }
                    }
                }
            }

            // Add the new vertices that appear in the skyline core
            for (int i = 0; i < coresAdded.size(); i++) {
                auto core = coresAdded[i];
                if (updatedMapCoCore.find(core) == updatedMapCoCore.end()) {
                    biSkylineCore tempSkylineCore = biSkylineCore();
                    tempSkylineCore.a = core.first;
                    tempSkylineCore.b = core.second;
                    updatedSkylineCores.push_back(tempSkylineCore);
                    updatedMapCoCore[core] = updatedSkylineCores.size() - 1;
                    originTree.Insert(updatedSkylineCores.size() - 1,
                                      bounds(tempSkylineCore.a, tempSkylineCore.b, 0, 0));
                }
                int coreId = updatedMapCoCore[core];
                updatedSkylineCores[coreId].vertices.push_back(vertex.vertexID);
                updatedSkylineCores[coreId].addedVertices.push_back(vertex.vertexID);
            }
        }
        end_time = clock();
        total_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
        map<int, vector<pair<int, int>>> mapCoreNode;  // coreID -> nodes
        updatedResortCores(originTree, updatedMapCoNode, updatedNodeMatrix, updatedMapCoCore, updatedCoreMatrix,
                           updatedSkylineCores, mapCoreNode);

        start_time = clock();
        set<int> changedCoreId;
        set<pair<int, int>> changedNodeCo;
        //// Update the connectivity of the vertices in the core and Rtree Node
        for (int tempA = alphaMax; tempA > 0; tempA--) {
            if (tempA < alphaMax && alphaBetaMax[tempA] == alphaBetaMax[tempA + 1]) {
                bool isBlank = true;
                for (int i = 1; i <= alphaBetaMax[tempA]; i++) {
                    if (updatedCoreMatrix[tempA][i] != -1 || updatedNodeMatrix[tempA][i] != -1) {
                        isBlank = false;
                        break;
                    }
                }
                if (isBlank) {
                    continue;
                }
            }
            UnionFind UF(uNum + vNum);
            vector<bool> nodes(uNum + vNum);
            std::fill_n(nodes.begin(), nodes.size(), false);
            BoundingBox bound;
            // check if there is new added core
            bool start = false;
            set<int> oldVertices;
            set<int> newVertices;
            for (int tempB = alphaBetaMax[tempA]; tempB > 0; tempB--) {
                vector<int> oldCores;
                vector<int> newCores;
                vector<int> tempResults;
                BoundingBox tempBound = bounds(tempA, tempB, alphaMax, 0);
                oldTree.QueryCore(RTree::AcceptEnclosing(tempBound), Visitor(), oldCores);
                newTree.QueryCore(RTree::AcceptEnclosing(tempBound), Visitor(), newCores);
                if (newCores.size() == 0)
                    continue;
                for (auto core: oldCores) {
                    for (auto vertex: oldSkylineCore[core].vertices) {
                        oldVertices.insert(vertex);
                    }
                }
                for (auto core: newCores) {
                    for (auto vertex: newSkylineCore[core].vertices) {
                        newVertices.insert(vertex);
                    }
                }

                std::set<int> addedVertices;
                std::set_difference(newVertices.begin(), newVertices.end(), oldVertices.begin(), oldVertices.end(),
                                    std::inserter(addedVertices, addedVertices.begin()));
                if (addedVertices.size() != 0) {
                    //// reconstruct the subgraphs
                    vector<int> cores;
                    if (!start) {
                        bound = bounds(tempA, tempB, alphaMax, betaMax);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        //// hash table
                        for (int tempA2 = tempA; tempA2 <= alphaMax; tempA2++) {
                            if (tempB > alphaBetaMax[tempA2]) {
                                break;
                            }
                            int coreId = updatedCoreMatrix[tempA2][tempB];
                            if (coreId != -1) {
                                cores.push_back(coreId);
                            }
                        }
                    }
                    vector<pair<int, int>> tempEdges;
                    unordered_map<int, int> newVertexCore;
                    for (int j = 0; j < cores.size(); j++) {
                        int tempCoreId = cores[j];
                        for (int m = 0; m < updatedSkylineCores[tempCoreId].subGraphs.size(); m++) {
                            vector<int> &vertices = updatedSkylineCores[tempCoreId].subGraphs[m];
                            if (vertices.size() > 0) {
                                auto it = vertices.begin();
                                int firstVertex = *it;
                                if (!nodes[firstVertex]) {
                                    nodes[firstVertex] = true;
                                }
                                newVertexCore[*it] = tempCoreId;
                                it++;
                                while (it != vertices.end()) {
                                    UF.Unite(firstVertex, *it);
                                    if (!nodes[*it]) {
                                        nodes[*it] = true;
                                        newVertexCore[*it] = tempCoreId;
                                    }
                                    ++it;
                                }
                            }
                        }
                        for (auto &edge: updatedSkylineCores[tempCoreId].connectEdges) {
                            tempEdges.push_back(edge);
                        }
                    }
                    for (auto &edge: tempEdges) {
                        UF.Unite(edge.first, edge.second);
                    }
//                    auto it = addedVertices.begin();
//                    nodes[*it] = true;
//                    int firstVertex = *it;
//                    it++;
//                    while (it != addedVertices.end()) {
//                        UF.Unite(firstVertex, *it);
//                        nodes[*it] = true;
//                        it++;
//                    }
                    // Check whether the new point can change the connectivity
                    vector<pair<int, int>> extraEdges;
                    auto it = addedVertices.begin();
                    while (it != addedVertices.end()) {
                        int vertexId = *it;
                        int root1 = UF.Find(vertexId);
                        if (vertexId < uNum) {
                            for (int node : g.neighbor_v1[vertexId]) {
                                node += uNum;
                                if (nodes[node]) {
                                    int root2 = UF.Find(node);
                                    if (root1 != root2 && addedVertices.find(node) == addedVertices.end()) {
                                        extraEdges.push_back({vertexId, node});
                                    }
                                    UF.Unite(vertexId, node);
                                }
                            }
                        } else {
                            for (int node : g.neighbor_v2[vertexId - uNum]) {
                                if (nodes[node]) {
                                    int root2 = UF.Find(node);
                                    if (root1 != root2 && addedVertices.find(node) == addedVertices.end()) {
                                        extraEdges.push_back({vertexId, node});
                                    }
                                    UF.Unite(vertexId, node);
                                }
                            }
                        }
                        it++;
                    }

                    if (updatedCoreMatrix[tempA][tempB] != -1) {
                        // updated the current core
                        bool changed = false;
                        int tempCoreID = updatedCoreMatrix[tempA][tempB];
                        map<int, vector<int>> subGraphRoot;  // root -> subgraph id
                        for (int j = 0; j < updatedSkylineCores[tempCoreID].subGraphs.size(); j++) {
                            int tempRoot = UF.Find(updatedSkylineCores[tempCoreID].subGraphs[j][0]);
                            if (subGraphRoot.find(tempRoot) != subGraphRoot.end()) {
                                changed = true;
                            }
                            subGraphRoot[UF.Find(updatedSkylineCores[tempCoreID].subGraphs[j][0])].push_back(j);
                        }
                        map<int, int> newSubGraphRoot;
                        if (changed) {
                            vector<vector<int>> tempSubGraphs;
                            int newGraphId = 0;
                            for (auto graphs: subGraphRoot) {
                                vector<int> tempSubgraph;
                                for (auto graphID: graphs.second) {
                                    tempSubgraph.insert(tempSubgraph.end(),
                                                        updatedSkylineCores[tempCoreID].subGraphs[graphID].begin(),
                                                        updatedSkylineCores[tempCoreID].subGraphs[graphID].end());
                                }
                                newSubGraphRoot[UF.Find(tempSubgraph[0])] = newGraphId;
                                newGraphId++;
                                tempSubGraphs.push_back(tempSubgraph);
                            }
                            updatedSkylineCores[tempCoreID].subGraphs = tempSubGraphs;
                            if (updatedSkylineCores[tempCoreID].addedVertices.size() != 0) {
                                for (auto tempVertex: updatedSkylineCores[tempCoreID].addedVertices) {
                                    if (newSubGraphRoot.find(UF.Find(tempVertex)) == newSubGraphRoot.end()) {
                                        newSubGraphRoot[UF.Find(
                                                tempVertex)] = updatedSkylineCores[tempCoreID].subGraphs.size();
                                        vector<int> tempGraph;
                                        tempGraph.push_back(tempVertex);
                                        updatedSkylineCores[tempCoreID].subGraphs.push_back(tempGraph);
                                    } else {
                                        updatedSkylineCores[tempCoreID].subGraphs[newSubGraphRoot[UF.Find(
                                                tempVertex)]].push_back(tempVertex);
                                    }
                                    newVertexCore[tempVertex] = tempCoreID;
                                }
                            }
//                            changedCoreId.insert(tempCoreID);
                            for (auto nodeCore: mapCoreNode[tempCoreID]) {
                                changedNodeCo.insert(nodeCore);
                            }
                        }
                    }
                    // update tempEdges
                    if (extraEdges.size() != 0) {
                        for (auto &tempEdge: extraEdges) {
                            int tempCoreId = newVertexCore[tempEdge.first];
                            updatedSkylineCores[tempCoreId].connectEdges.push_back(tempEdge);
                        }
                    }
                } else {
                    //// reconstruct the subgraphs
                    vector<int> cores;
                    if (!start) {
                        bound = bounds(tempA, tempB, alphaMax, betaMax);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        //// hash table
                        for (int tempA2 = tempA; tempA2 <= alphaMax; tempA2++) {
                            if (tempB > alphaBetaMax[tempA2]) {
                                break;
                            }
                            int coreId = updatedCoreMatrix[tempA2][tempB];
                            if (coreId != -1) {
                                cores.push_back(coreId);
                            }
                        }
                    }
                    vector<pair<int, int>> tempEdges;
                    unordered_map<int, int> newVertexCore;
                    for (int j = 0; j < cores.size(); j++) {
                        int tempCoreId = cores[j];
                        for (int m = 0; m < updatedSkylineCores[tempCoreId].subGraphs.size(); m++) {
                            vector<int> &vertices = updatedSkylineCores[tempCoreId].subGraphs[m];
                            if (vertices.size() > 0) {
                                auto it = vertices.begin();
                                int firstVertex = *it;
                                if (!nodes[firstVertex]) {
                                    nodes[firstVertex] = true;
                                }
                                newVertexCore[*it] = tempCoreId;
                                it++;
                                while (it != vertices.end()) {
                                    UF.Unite(firstVertex, *it);
                                    if (!nodes[*it]) {
                                        nodes[*it] = true;
                                        newVertexCore[*it] = tempCoreId;
                                    }
                                    ++it;
                                }
                            }
                        }
                        for (auto edge: updatedSkylineCores[tempCoreId].connectEdges) {
                            tempEdges.push_back(edge);
                        }
                    }
                    for (auto &edge: tempEdges) {
                        UF.Unite(edge.first, edge.second);
                    }
//                    cout << tempCase.fromVertex << " " << tempCase.toVertex << " ";
                    bool change = false;
                    for (auto &vertexPair: tempCase.vertexPairs) {
                        if (UF.Find(vertexPair.first) != UF.Find(vertexPair.second)) {
                            change = true;
                            UF.Unite(vertexPair.first, vertexPair.second);
                            updatedSkylineCores[newVertexCore[vertexPair.first]].connectEdges.push_back(vertexPair);
                            for (auto nodeCore: mapCoreNode[newVertexCore[vertexPair.first]]) {
                                changedNodeCo.insert(nodeCore);
                            }
                        }
                    }
                    if (!change && newVertices.size() == updatedVertices.size()) {
                        break;
                    }
                }

                // updated R-tree node
                if (changedNodeCo.find({tempA, tempB}) != changedNodeCo.end()) {
                    auto RNodes = updatedMapCoNode.find({tempA, tempB})->second;
                    vector<bool> tempNodes;
                    tempNodes.resize(uNum + vNum);
                    vector<int> subgraphRoot;
                    vector<int> tempTag;
                    for (auto RNode: RNodes) {
                        int coreLeftId = RNode->left;
                        int coreRightId = RNode->right;
                        // get vertices in the range of R-Node and get the tag
                        for (int j = coreLeftId; j <= coreRightId; j++) {
                            for (auto &tempSubgraph: updatedSkylineCores[j].subGraphs) {
                                for (auto &tempVertex: tempSubgraph) {
                                    tempNodes[tempVertex] = true;
                                }
                                int tag = -1;
                                int tempRoot = UF.Find(tempSubgraph[0]);
                                for (int k = 0; k < subgraphRoot.size(); k++) {
                                    if (subgraphRoot[k] == tempRoot) {
                                        tag = k;
                                        break;
                                    }
                                }
                                if (tag == -1) {
                                    subgraphRoot.push_back(tempRoot);
                                    tag = subgraphRoot.size() - 1;
                                }
                                tempTag.push_back(tag);
                            }
                        }
                        RNode->tag = tempTag;
                        // get out edges of R-Node
                        RNode->outEdges.clear();
                        for (int j = coreLeftId; j <= coreRightId; j++) {
                            for (auto &tempEdge: updatedSkylineCores[j].connectEdges) {
                                if (!tempNodes[tempEdge.first] || !tempNodes[tempEdge.second]) {
                                    RNode->outEdges.push_back(tempEdge);
                                }
                            }
                        }
                    }
                    changedNodeCo.erase({tempA, tempB});
                }
            }
        }
        end_time = clock();
//        cout << changedCoreId.size() << endl;
        total_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
    }
    cout << "Update time : " << total_time / (double) testCase.size() << "s" << endl;
//    cout << "Update time : " << total_time / (double) valid << endl;
}

void addOneMatrix(RTree tree, BiGraph &g, int addNum) {
    //// load the file
    string addOnePath = "./ab-core/" + dataset + "/" + dataset + "-add";
    addOnePath += to_string(addNum);
    ifstream addOneFile(addOnePath);
    string line;
    vector<addOneEdge> testCase;
    while (std::getline(addOneFile, line)) {
        if (line.size() <= 0)
            break;
        if (line[0] == 'A') {
            int pos = line.find(":");
            line = line.substr(pos + 2);
            istringstream iss(line);
            int fromVertex, toVertex;
            addOneEdge tempAdd;
            while (iss >> fromVertex >> toVertex) {
//                cout << fromVertex << " " << toVertex << endl;
                tempAdd.vertexPairs.push_back({fromVertex, toVertex});
            }
            testCase.push_back(tempAdd);
        } else {
            string delimiter1 = ":";
            string delimiter2 = "|";
            int pos1 = line.find(delimiter1);
            int pos2 = line.find(delimiter2);
            string strID = line.substr(0, pos1);
            string strCoordinatesNew = line.substr(pos1 + 1, pos2 - pos1 - 1);
            string strCoordinatesOld = line.substr(pos2 + 1);

            UpdatedVertex tempVertex;
            tempVertex.vertexID = stoi(strID);

            string delimiter = ")(";

            //// Process the new coordinates
            strCoordinatesNew = strCoordinatesNew.substr(1, strCoordinatesNew.length() - 2);
            size_t pos = 0;
            std::string token;
            while ((pos = strCoordinatesNew.find(delimiter)) != std::string::npos) {
                token = strCoordinatesNew.substr(0, pos); // Extract a coordinate pair string
                // Extract the x and y strings
                size_t commaPos = token.find(",");
                std::string x_str = token.substr(0, commaPos);
                std::string y_str = token.substr(commaPos + 1);
                // Convert strings to integers
                int x = std::stoi(x_str);
                int y = std::stoi(y_str);
                tempVertex.newCores.push_back({x, y});
                strCoordinatesNew.erase(0, pos + delimiter.length()); // Remove the processed part
            }
            int commaPos = strCoordinatesNew.find(",");
            std::string x_str = strCoordinatesNew.substr(0, commaPos);
            std::string y_str = strCoordinatesNew.substr(commaPos + 1);
            int x = std::stoi(x_str);
            int y = std::stoi(y_str);
            tempVertex.newCores.push_back({x, y});

            //// Process the old coordinates
            strCoordinatesOld = strCoordinatesOld.substr(1, strCoordinatesOld.length() - 2);
            pos = 0;
            while ((pos = strCoordinatesOld.find(delimiter)) != std::string::npos) {
                token = strCoordinatesOld.substr(0, pos); // Extract a coordinate pair string
                // Extract the x and y strings
                commaPos = token.find(",");
                x_str = token.substr(0, commaPos);
                y_str = token.substr(commaPos + 1);
                // Convert strings to integers
                x = std::stoi(x_str);
                y = std::stoi(y_str);
                tempVertex.oldCores.push_back({x, y});
                strCoordinatesOld.erase(0, pos + delimiter.length()); // Remove the processed part
            }
            commaPos = strCoordinatesOld.find(",");
            x_str = strCoordinatesOld.substr(0, commaPos);
            y_str = strCoordinatesOld.substr(commaPos + 1);
            x = std::stoi(x_str);
            y = std::stoi(y_str);
            tempVertex.oldCores.push_back({x, y});
            testCase[testCase.size() - 1].updatedVertices.push_back(tempVertex);
        }
    }
    addOneFile.close();

    int cnt = 1;
    clock_t start_time, end_time;
    double total_time = 0;
    for (auto tempCase: testCase) {
//        cout << cnt << "/" << testCase.size() << endl;
        cnt++;
        // Get skyline cores in updated vertices
        auto &updatedVertices = tempCase.updatedVertices;
        if (updatedVertices.size() == 0)
            continue;
        //// build temp index for test
        auto updatedSkylineCores = skylineCores;
        map<pair<int, int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>> updatedMapCoNode;
        map<pair<int, int>, int> updatedMapCoCore = mapCoCore;
        vector<vector<int>> updatedCoreMatrix = coreMatrix;
        vector<vector<int>> updatedNodeMatrix = nodeMatrix;
        RTree originTree;
        for (int i = 0; i < updatedSkylineCores.size(); i++) {
            int xx = updatedSkylineCores[i].a;
            int yy = updatedSkylineCores[i].b;
            originTree.Insert(i, bounds(xx, yy, 0, 0));
        }

        map<pair<int, int>, vector<int>> mapOldSkylineCore;
        map<pair<int, int>, vector<int>> mapNewSkylineCore;
        for (auto &tempVertex: updatedVertices) {
            for (int i = 0; i < tempVertex.oldCores.size(); i++) {
                mapOldSkylineCore[{tempVertex.oldCores[i].first, tempVertex.oldCores[i].second}].push_back(
                        tempVertex.vertexID);
            }
            for (int i = 0; i < tempVertex.newCores.size(); i++) {
                mapNewSkylineCore[{tempVertex.newCores[i].first, tempVertex.newCores[i].second}].push_back(
                        tempVertex.vertexID);
            }
        }
        // build 2 matrix for updated vertices
        vector<vector<int>> oldMatrix;
        vector<vector<int>> newMatrix;
        oldMatrix.resize(alphaMax + 1);
        newMatrix.resize(alphaMax + 1);
        int updatedAlphaMax = 0;
        vector<int> updatedKVertexNum(alphaMax + 1);
        for (int i = 0; i < alphaMax + 1; i++) {
            oldMatrix[i].resize(alphaBetaMax[i] + 1);
            newMatrix[i].resize(alphaBetaMax[i] + 1);
            fill_n(oldMatrix[i].begin(), oldMatrix[i].size(), -1);
            fill_n(newMatrix[i].begin(), newMatrix[i].size(), -1);
        }
        vector<biSkylineCore> oldSkylineCore;
        vector<biSkylineCore> newSkylineCore;
        map<pair<int, int>, int> mapOldCoreId;
        map<pair<int, int>, int> mapNewCoreId;
        int index1 = 0;
        int index2 = 0;
        for (auto &core: mapOldSkylineCore) {
            int xx = core.first.first;  // k
            int yy = core.first.second; // l
            biSkylineCore tempSkylineCore = biSkylineCore();
            tempSkylineCore.a = core.first.first;
            tempSkylineCore.b = core.first.second;
            tempSkylineCore.vertices = core.second;
            std::sort(tempSkylineCore.vertices.begin(), tempSkylineCore.vertices.end());
            oldSkylineCore.push_back(tempSkylineCore);
            mapOldCoreId[{tempSkylineCore.a, tempSkylineCore.b}] = index1;
            oldMatrix[xx][yy] = index1;
            index1++;
        }
        for (auto &core: mapNewSkylineCore) {
            int xx = core.first.first;  // k
            int yy = core.first.second; // l
            biSkylineCore tempSkylineCore = biSkylineCore();
            tempSkylineCore.a = core.first.first;
            tempSkylineCore.b = core.first.second;
            tempSkylineCore.vertices = core.second;
            std::sort(tempSkylineCore.vertices.begin(), tempSkylineCore.vertices.end());
            newSkylineCore.push_back(tempSkylineCore);
            mapNewCoreId[{tempSkylineCore.a, tempSkylineCore.b}] = index2;
            newMatrix[xx][yy] = index2;
            if (oldMatrix[xx][yy] != -1) {
                if (oldSkylineCore[oldMatrix[xx][yy]].vertices != tempSkylineCore.vertices) {
                    updatedAlphaMax = max(updatedAlphaMax, xx);
                }
            } else {
                updatedAlphaMax = max(updatedAlphaMax, xx);
            }
            updatedAlphaMax = max(updatedAlphaMax, xx);
            index2++;
        }
//        cout << updatedAlphaMax << endl;
        int cnt = 0;
        set<int> tempUpdatedVertex;
        for (int tempA = alphaMax; tempA > 0; tempA--) {
            for (int tempB = alphaBetaMax[tempA]; tempB > 0; tempB--) {
                if (oldMatrix[tempA][tempB] != -1) {
                    for (auto &vertex : oldSkylineCore[oldMatrix[tempA][tempB]].vertices) {
                        tempUpdatedVertex.insert(vertex);
                    }
                }
            }
            updatedKVertexNum[tempA] = tempUpdatedVertex.size();
        }

        start_time = clock();
        // updated original skyline core
        for (auto vertex: updatedVertices) {
            // find the changed skyline cores of vertex
            std::sort(vertex.newCores.begin(), vertex.newCores.end());
            std::sort(vertex.oldCores.begin(), vertex.oldCores.end());
            vector<pair<int, int>> coresAdded;
            vector<pair<int, int>> coresDeleted;
            std::set_difference(vertex.newCores.begin(), vertex.newCores.end(),
                                vertex.oldCores.begin(), vertex.oldCores.end(),
                                std::back_inserter(coresAdded));
            std::set_difference(vertex.oldCores.begin(), vertex.oldCores.end(),
                                vertex.newCores.begin(), vertex.newCores.end(),
                                std::back_inserter(coresDeleted));
            // Delete vertices in skyline cores that no longer exist
            for (int i = 0; i < coresDeleted.size(); i++) {
                int coreId = updatedMapCoCore[coresDeleted[i]];
                // remove the vertex in the old core
                auto oldEnd = std::remove(updatedSkylineCores[coreId].vertices.begin(),
                                          updatedSkylineCores[coreId].vertices.end(),
                                          vertex.vertexID);
                updatedSkylineCores[coreId].vertices.erase(oldEnd, updatedSkylineCores[coreId].vertices.end());
                if (updatedSkylineCores[coreId].vertices.size() == 0) {
                    BoundingBox bound = bounds(coresDeleted[i].first, coresDeleted[i].second, 0, 0);
                    originTree.RemoveBoundedArea(bound);
                    updatedMapCoCore.erase({coresDeleted[i].first, coresDeleted[i].second});
                    updatedCoreMatrix[coresDeleted[i].first][coresDeleted[i].second] = -1;
                }
                for (int j = 0; j < updatedSkylineCores[coreId].subGraphs.size(); j++) {
                    auto iter = std::find(updatedSkylineCores[coreId].subGraphs[j].begin(),
                                          updatedSkylineCores[coreId].subGraphs[j].end(), vertex.vertexID);
                    if (iter != updatedSkylineCores[coreId].subGraphs[j].end()) {
                        updatedSkylineCores[coreId].subGraphs[j].erase(iter);
                        if (updatedSkylineCores[coreId].subGraphs[j].size() == 0) {
                            updatedSkylineCores[coreId].subGraphs.erase(
                                    updatedSkylineCores[coreId].subGraphs.begin() + j);
                        }
                    }
                }
            }

            // Add the new vertices that appear in the skyline core
            for (int i = 0; i < coresAdded.size(); i++) {
                auto core = coresAdded[i];
                if (updatedMapCoCore.find(core) == updatedMapCoCore.end()) {
                    biSkylineCore tempSkylineCore = biSkylineCore();
                    tempSkylineCore.a = core.first;
                    tempSkylineCore.b = core.second;
                    updatedSkylineCores.push_back(tempSkylineCore);
                    updatedMapCoCore[core] = updatedSkylineCores.size() - 1;
                    originTree.Insert(updatedSkylineCores.size() - 1,
                                      bounds(tempSkylineCore.a, tempSkylineCore.b, 0, 0));
                }
                int coreId = updatedMapCoCore[core];
                updatedSkylineCores[coreId].vertices.push_back(vertex.vertexID);
                updatedSkylineCores[coreId].addedVertices.push_back(vertex.vertexID);
            }
        }
        end_time = clock();
        total_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
        map<int, vector<pair<int, int>>> mapCoreNode;  // coreID -> nodes
        updatedResortCores(originTree, updatedMapCoNode, updatedNodeMatrix, updatedMapCoCore, updatedCoreMatrix,
                           updatedSkylineCores, mapCoreNode);

        start_time = clock();
        set<int> changedCoreId;
        set<pair<int, int>> changedNodeCo;

        UnionFind UF(uNum + vNum);
        vector<bool> nodes(uNum + vNum);
        std::fill_n(nodes.begin(), nodes.size(), false);
        vector<int> preResult;
        vector<bool> oldFlag(uNum + vNum);
        vector<int> preOld;
        cout << "updated Alpha Max:" << updatedAlphaMax << endl;
        //// Update the connectivity of the vertices in the core and Rtree Node
        for (int tempK = updatedAlphaMax; tempK > 0; tempK--) {
            if (tempK < alphaMax && alphaBetaMax[tempK] == alphaBetaMax[tempK + 1]) {
                bool isBlank = true;
                for (int i = 1; i <= alphaBetaMax[tempK]; i++) {
                    if (updatedCoreMatrix[tempK][i] != -1 || updatedNodeMatrix[tempK][i] != -1) {
                        isBlank = false;
                        break;
                    }
                }
                if (isBlank) {
                    continue;
                }
            }
            UF.Init(preResult);
            for (auto &element: preResult) {
                nodes[element] = false;
            }
            preResult.clear();
            BoundingBox bound;
            // check if there is new added core
            bool start = false;
//            set<int> oldVertices;
//            set<int> newVertices;
            int oldCnt = 0;
            set<int> addedVertices;
            for (auto &element : preOld) {
                oldFlag[element] = false;
            }
            preOld.clear();
            for (int tempB = alphaBetaMax[tempK]; tempB > 0; tempB--) {
                vector<int> oldCores;
                vector<int> newCores;
                vector<int> tempResults;
                for (int tempA2 = tempK; tempA2 <= alphaMax; tempA2++) {
                    if (tempB > alphaBetaMax[tempA2]) {
                        break;
                    }
                    int coreId = oldMatrix[tempA2][tempB];
                    if (coreId != -1) {
                        oldCores.push_back(coreId);
                    }
                    coreId = newMatrix[tempA2][tempB];
                    if (coreId != -1) {
                        newCores.push_back(coreId);
                    }
                }
                if (newCores.size() == 0)
                    continue;
                for (auto core: oldCores) {
                    for (auto vertex: oldSkylineCore[core].vertices) {
                        if (!oldFlag[vertex]) {
                            oldFlag[vertex] = true;
                            preOld.push_back(vertex);
                            addedVertices.erase(vertex);
                            oldCnt++;
                        }
                    }
                }
                for (auto core: newCores) {
                    for (auto vertex: newSkylineCore[core].vertices) {
                        if (!oldFlag[vertex]) {
                            addedVertices.insert(vertex);
                        }
                    }
                }

//                std::set<int> addedVertices;
//                std::set_difference(newVertices.begin(), newVertices.end(), oldVertices.begin(), oldVertices.end(),
//                                    std::inserter(addedVertices, addedVertices.begin()));
                if (addedVertices.size() != 0) {
                    //// reconstruct the subgraphs
                    vector<int> cores;
                    if (!start) {
                        bound = bounds(tempK, tempB, alphaMax, betaMax);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        //// hash table
                        for (int tempA2 = tempK; tempA2 <= alphaMax; tempA2++) {
                            if (tempB > alphaBetaMax[tempA2]) {
                                break;
                            }
                            int coreId = updatedCoreMatrix[tempA2][tempB];
                            if (coreId != -1) {
                                cores.push_back(coreId);
                            }
                        }
                    }
                    vector<pair<int, int>> tempEdges;
                    unordered_map<int, int> newVertexCore;
                    for (int j = 0; j < cores.size(); j++) {
                        int tempCoreId = cores[j];
                        for (int m = 0; m < updatedSkylineCores[tempCoreId].subGraphs.size(); m++) {
                            vector<int> &vertices = updatedSkylineCores[tempCoreId].subGraphs[m];
                            if (vertices.size() > 0) {
                                auto it = vertices.begin();
                                int firstVertex = *it;
                                if (!nodes[firstVertex]) {
                                    nodes[firstVertex] = true;
                                }
                                newVertexCore[*it] = tempCoreId;
                                it++;
                                while (it != vertices.end()) {
                                    UF.Unite(firstVertex, *it);
                                    if (!nodes[*it]) {
                                        nodes[*it] = true;
                                        preResult.push_back(*it);
                                        newVertexCore[*it] = tempCoreId;
                                    }
                                    ++it;
                                }
                            }
                        }
                        for (auto &edge: updatedSkylineCores[tempCoreId].connectEdges) {
                            tempEdges.push_back(edge);
                        }
                    }
                    for (auto &edge: tempEdges) {
                        UF.Unite(edge.first, edge.second);
                    }
//                    auto it = addedVertices.begin();
//                    nodes[*it] = true;
//                    int firstVertex = *it;
//                    it++;
//                    while (it != addedVertices.end()) {
//                        UF.Unite(firstVertex, *it);
//                        nodes[*it] = true;
//                        it++;
//                    }
                    // Check whether the new point can change the connectivity
                    vector<pair<int, int>> extraEdges;
                    auto it = addedVertices.begin();
                    while (it != addedVertices.end()) {
                        int vertexId = *it;
                        int root1 = UF.Find(vertexId);
                        if (vertexId < uNum) {
                            for (int node : g.neighbor_v1[vertexId]) {
                                node += uNum;
                                if (nodes[node]) {
                                    int root2 = UF.Find(node);
                                    if (root1 != root2 && addedVertices.find(node) == addedVertices.end()) {
                                        extraEdges.push_back({vertexId, node});
                                    }
                                    UF.Unite(vertexId, node);
                                }
                            }
                        } else {
                            for (int node : g.neighbor_v2[vertexId - uNum]) {
                                if (nodes[node]) {
                                    int root2 = UF.Find(node);
                                    if (root1 != root2 && addedVertices.find(node) == addedVertices.end()) {
                                        extraEdges.push_back({vertexId, node});
                                    }
                                    UF.Unite(vertexId, node);
                                }
                            }
                        }
                        it++;
                    }

                    if (updatedCoreMatrix[tempK][tempB] != -1) {
                        // updated the current core
                        bool changed = false;
                        int tempCoreID = updatedCoreMatrix[tempK][tempB];
                        map<int, vector<int>> subGraphRoot;  // root -> subgraph id
                        for (int j = 0; j < updatedSkylineCores[tempCoreID].subGraphs.size(); j++) {
                            int tempRoot = UF.Find(updatedSkylineCores[tempCoreID].subGraphs[j][0]);
                            if (subGraphRoot.find(tempRoot) != subGraphRoot.end()) {
                                changed = true;
                            }
                            subGraphRoot[UF.Find(updatedSkylineCores[tempCoreID].subGraphs[j][0])].push_back(j);
                        }
                        map<int, int> newSubGraphRoot;
                        if (changed) {
                            vector<vector<int>> tempSubGraphs;
                            int newGraphId = 0;
                            for (auto graphs: subGraphRoot) {
                                vector<int> tempSubgraph;
                                for (auto graphID: graphs.second) {
                                    tempSubgraph.insert(tempSubgraph.end(),
                                                        updatedSkylineCores[tempCoreID].subGraphs[graphID].begin(),
                                                        updatedSkylineCores[tempCoreID].subGraphs[graphID].end());
                                }
                                newSubGraphRoot[UF.Find(tempSubgraph[0])] = newGraphId;
                                newGraphId++;
                                tempSubGraphs.push_back(tempSubgraph);
                            }
                            updatedSkylineCores[tempCoreID].subGraphs = tempSubGraphs;
                            if (updatedSkylineCores[tempCoreID].addedVertices.size() != 0) {
                                for (auto tempVertex: updatedSkylineCores[tempCoreID].addedVertices) {
                                    if (newSubGraphRoot.find(UF.Find(tempVertex)) == newSubGraphRoot.end()) {
                                        newSubGraphRoot[UF.Find(
                                                tempVertex)] = updatedSkylineCores[tempCoreID].subGraphs.size();
                                        vector<int> tempGraph;
                                        tempGraph.push_back(tempVertex);
                                        updatedSkylineCores[tempCoreID].subGraphs.push_back(tempGraph);
                                    } else {
                                        updatedSkylineCores[tempCoreID].subGraphs[newSubGraphRoot[UF.Find(
                                                tempVertex)]].push_back(tempVertex);
                                    }
                                    newVertexCore[tempVertex] = tempCoreID;
                                }
                            }
//                            changedCoreId.insert(tempCoreID);
                            for (auto nodeCore: mapCoreNode[tempCoreID]) {
                                changedNodeCo.insert(nodeCore);
                            }
                        }
                    }
                    // update tempEdges
                    if (extraEdges.size() != 0) {
                        for (auto &tempEdge: extraEdges) {
                            int tempCoreId = newVertexCore[tempEdge.first];
                            updatedSkylineCores[tempCoreId].connectEdges.push_back(tempEdge);
                        }
                    }
                } else {
                    //// reconstruct the subgraphs
                    vector<int> cores;
                    if (!start) {
                        bound = bounds(tempK, tempB, alphaMax, betaMax);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        //// hash table
                        for (int tempA2 = tempK; tempA2 <= alphaMax; tempA2++) {
                            if (tempB > alphaBetaMax[tempA2]) {
                                break;
                            }
                            int coreId = updatedCoreMatrix[tempA2][tempB];
                            if (coreId != -1) {
                                cores.push_back(coreId);
                            }
                        }
                    }
                    vector<pair<int, int>> tempEdges;
                    unordered_map<int, int> newVertexCore;
                    for (int j = 0; j < cores.size(); j++) {
                        int tempCoreId = cores[j];
                        for (int m = 0; m < updatedSkylineCores[tempCoreId].subGraphs.size(); m++) {
                            vector<int> &vertices = updatedSkylineCores[tempCoreId].subGraphs[m];
                            if (vertices.size() > 0) {
                                auto it = vertices.begin();
                                int firstVertex = *it;
                                if (!nodes[firstVertex]) {
                                    nodes[firstVertex] = true;
                                }
                                newVertexCore[*it] = tempCoreId;
                                it++;
                                while (it != vertices.end()) {
                                    UF.Unite(firstVertex, *it);
                                    if (!nodes[*it]) {
                                        nodes[*it] = true;
                                        preResult.push_back(*it);
                                        newVertexCore[*it] = tempCoreId;
                                    }
                                    ++it;
                                }
                            }
                        }
                        for (auto edge: updatedSkylineCores[tempCoreId].connectEdges) {
                            tempEdges.push_back(edge);
                        }
                    }
                    for (auto &edge: tempEdges) {
                        UF.Unite(edge.first, edge.second);
                    }
//                    cout << tempCase.fromVertex << " " << tempCase.toVertex << " ";
                    bool change = false;
                    for (auto &vertexPair: tempCase.vertexPairs) {
                        if (UF.Find(vertexPair.first) != UF.Find(vertexPair.second)) {
                            change = true;
                            UF.Unite(vertexPair.first, vertexPair.second);
                            updatedSkylineCores[newVertexCore[vertexPair.first]].connectEdges.push_back(vertexPair);
                            for (auto nodeCore: mapCoreNode[newVertexCore[vertexPair.first]]) {
                                changedNodeCo.insert(nodeCore);
                            }
                        }
                    }
                    if (!change && oldCnt >= updatedKVertexNum[tempK]) {
                        break;
                    }
                }

                // updated R-tree node
                if (changedNodeCo.find({tempK, tempB}) != changedNodeCo.end()) {
                    auto RNodes = updatedMapCoNode.find({tempK, tempB})->second;
                    vector<bool> tempNodes;
                    tempNodes.resize(uNum + vNum);
                    vector<int> subgraphRoot;
                    vector<int> tempTag;
                    for (auto RNode: RNodes) {
                        int coreLeftId = RNode->left;
                        int coreRightId = RNode->right;
                        // get vertices in the range of R-Node and get the tag
                        for (int j = coreLeftId; j <= coreRightId; j++) {
                            for (auto &tempSubgraph: updatedSkylineCores[j].subGraphs) {
                                for (auto &tempVertex: tempSubgraph) {
                                    tempNodes[tempVertex] = true;
                                }
                                int tag = -1;
                                int tempRoot = UF.Find(tempSubgraph[0]);
                                for (int k = 0; k < subgraphRoot.size(); k++) {
                                    if (subgraphRoot[k] == tempRoot) {
                                        tag = k;
                                        break;
                                    }
                                }
                                if (tag == -1) {
                                    subgraphRoot.push_back(tempRoot);
                                    tag = subgraphRoot.size() - 1;
                                }
                                tempTag.push_back(tag);
                            }
                        }
                        RNode->tag = tempTag;
                        // get out edges of R-Node
                        RNode->outEdges.clear();
                        for (int j = coreLeftId; j <= coreRightId; j++) {
                            for (auto &tempEdge: updatedSkylineCores[j].connectEdges) {
                                if (!tempNodes[tempEdge.first] || !tempNodes[tempEdge.second]) {
                                    RNode->outEdges.push_back(tempEdge);
                                }
                            }
                        }
                    }
                    changedNodeCo.erase({tempK, tempB});
                }
            }
        }
        end_time = clock();
//        cout << changedCoreId.size() << endl;
        total_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
    }
    cout << "Update time : " << total_time / (double) testCase.size() << "s" << endl;
//    cout << "Update time : " << total_time / (double) valid << endl;
}

void deleteOne(RTree tree, BiGraph &g, int deleteNum) {
    //// load the file
    string deleteOnePath = "./ab-core/" + dataset + "/" + dataset + "-delete";
    deleteOnePath += to_string(deleteNum);
    ifstream deleteOneFile(deleteOnePath);
    string line;
    vector<deleteOneEdge> testCase;
    while (std::getline(deleteOneFile, line)) {
        if (line.size() <= 0)
            break;
        if (line[0] == 'D') {
            int pos = line.find(":");
            line = line.substr(pos + 2);
            istringstream iss(line);
            int fromVertex, toVertex;
            deleteOneEdge tempDelete;
            while (iss >> fromVertex >> toVertex) {
//                cout << fromVertex << " " << toVertex << endl;
                tempDelete.vertexPairs.push_back({fromVertex, toVertex});
            }
//            cout << tempDelete.vertexPairs.size() << endl;
            testCase.push_back(tempDelete);
        } else {
            string delimiter1 = ":";
            string delimiter2 = "|";
            int pos1 = line.find(delimiter1);
            int pos2 = line.find(delimiter2);
            string strID = line.substr(0, pos1);
            string strCoordinatesNew = line.substr(pos1 + 1, pos2 - pos1 - 1);
            string strCoordinatesOld = line.substr(pos2 + 1);

            UpdatedVertex tempVertex;
            tempVertex.vertexID = stoi(strID);

            string delimiter = ")(";

            //// Process the new coordinates
            if (strCoordinatesNew.size() > 0) {
                strCoordinatesNew = strCoordinatesNew.substr(1, strCoordinatesNew.length() - 2);
                size_t pos = 0;
                std::string token;
                while ((pos = strCoordinatesNew.find(delimiter)) != std::string::npos) {
                    token = strCoordinatesNew.substr(0, pos); // Extract a coordinate pair string
                    // Extract the x and y strings
                    size_t commaPos = token.find(",");
                    std::string x_str = token.substr(0, commaPos);
                    std::string y_str = token.substr(commaPos + 1);
                    // Convert strings to integers
                    int x = std::stoi(x_str);
                    int y = std::stoi(y_str);
                    tempVertex.newCores.push_back({x, y});
                    strCoordinatesNew.erase(0, pos + delimiter.length()); // Remove the processed part
                }
                int commaPos = strCoordinatesNew.find(",");
                std::string x_str = strCoordinatesNew.substr(0, commaPos);
                std::string y_str = strCoordinatesNew.substr(commaPos + 1);
                int x = std::stoi(x_str);
                int y = std::stoi(y_str);
                tempVertex.newCores.push_back({x, y});
            }

            //// Process the old coordinates
            if (strCoordinatesOld.size() > 0) {
                strCoordinatesOld = strCoordinatesOld.substr(1, strCoordinatesOld.length() - 2);
                size_t pos = 0;
                std::string token;
                string x_str;
                string y_str;
                int x;
                int y;
                size_t commaPos;
                while ((pos = strCoordinatesOld.find(delimiter)) != std::string::npos) {
                    token = strCoordinatesOld.substr(0, pos); // Extract a coordinate pair string
                    // Extract the x and y strings
                    commaPos = token.find(",");
                    x_str = token.substr(0, commaPos);
                    y_str = token.substr(commaPos + 1);
                    // Convert strings to integers
                    x = std::stoi(x_str);
                    y = std::stoi(y_str);
                    tempVertex.oldCores.push_back({x, y});
                    strCoordinatesOld.erase(0, pos + delimiter.length()); // Remove the processed part
                }
                commaPos = strCoordinatesOld.find(",");
                x_str = strCoordinatesOld.substr(0, commaPos);
                y_str = strCoordinatesOld.substr(commaPos + 1);
                x = std::stoi(x_str);
                y = std::stoi(y_str);
                tempVertex.oldCores.push_back({x, y});
                testCase[testCase.size() - 1].updatedVertices.push_back(tempVertex);
            }
        }
    }
    deleteOneFile.close();

    int cnt = 1;
    int valid = 0;
    clock_t start_time, end_time;
    double total_time = 0;
    for (auto tempCase: testCase) {
        double this_time = 0;
        cout << cnt << "/" << testCase.size() << endl;
        cnt++;
        // Get skyline cores in updated vertices
        auto &updatedVertices = tempCase.updatedVertices;
        if (updatedVertices.size() == 0)
            continue;
        valid++;
        //// build temp index for test
        auto updatedV1Graph = g.neighbor_v1;
        auto updatedV2Graph = g.neighbor_v2;
        for (auto &vertexPair: tempCase.vertexPairs) {
            updatedV1Graph[vertexPair.first].erase(
                    std::remove(updatedV1Graph[vertexPair.first].begin(),
                                updatedV1Graph[vertexPair.first].end(),
                                vertexPair.second), updatedV1Graph[vertexPair.first].end());
            updatedV2Graph[vertexPair.second].erase(
                    std::remove(updatedV2Graph[vertexPair.second].begin(), updatedV2Graph[vertexPair.second].end(),
                                vertexPair.first), updatedV2Graph[vertexPair.second].end());
        }

        auto updatedSkylineCores = skylineCores;
        map<pair<int, int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>> updatedMapCoNode;
        map<pair<int, int>, int> updatedMapCoCore = mapCoCore;
        vector<vector<int>> updatedCoreMatrix = coreMatrix;
        vector<vector<int>> updatedNodeMatrix = nodeMatrix;
        RTree originTree;
        for (int i = 0; i < updatedSkylineCores.size(); i++) {
            int xx = updatedSkylineCores[i].a;
            int yy = updatedSkylineCores[i].b;
            originTree.Insert(i, bounds(xx, yy, 0, 0));
        }

        map<pair<int, int>, vector<int>> mapOldSkylineCore;
        map<pair<int, int>, vector<int>> mapNewSkylineCore;
        for (auto &tempVertex: updatedVertices) {
            for (int i = 0; i < tempVertex.oldCores.size(); i++) {
                mapOldSkylineCore[{tempVertex.oldCores[i].first, tempVertex.oldCores[i].second}].push_back(
                        tempVertex.vertexID);
            }
            for (int i = 0; i < tempVertex.newCores.size(); i++) {
                mapNewSkylineCore[{tempVertex.newCores[i].first, tempVertex.newCores[i].second}].push_back(
                        tempVertex.vertexID);
            }
        }
        start_time = clock();
        // build 2 R-tree for updated vertices
        vector<biSkylineCore> oldSkylineCore;
        vector<biSkylineCore> newSkylineCore;
        RTree oldTree;
        RTree newTree;
        map<pair<int, int>, int> mapOldCoreId;
        map<pair<int, int>, int> mapNewCoreId;
        int index1 = 0;
        int index2 = 0;
        for (auto &core: mapOldSkylineCore) {
            int xx = core.first.first;  // k
            int yy = core.first.second; // l
            oldTree.Insert(index1, bounds(xx, yy, 0, 0));
            biSkylineCore tempSkylineCore = biSkylineCore();
            tempSkylineCore.a = core.first.first;
            tempSkylineCore.b = core.first.second;
            tempSkylineCore.vertices = core.second;
            oldSkylineCore.push_back(tempSkylineCore);
            mapOldCoreId[{tempSkylineCore.a, tempSkylineCore.b}] = index1;
            index1++;
        }
        for (auto &core: mapNewSkylineCore) {
            int xx = core.first.first;  // k
            int yy = core.first.second; // l
            newTree.Insert(index2, bounds(xx, yy, 0, 0));
            biSkylineCore tempSkylineCore = biSkylineCore();
            tempSkylineCore.a = core.first.first;
            tempSkylineCore.b = core.first.second;
            tempSkylineCore.vertices = core.second;
            newSkylineCore.push_back(tempSkylineCore);
            mapNewCoreId[{tempSkylineCore.a, tempSkylineCore.b}] = index2;
            index2++;
        }
        long long size = 0;
        oldTree.indexSize(RTree::AcceptAny(), Visitor(), size, mapOldCoreId);
        newTree.indexSize(RTree::AcceptAny(), Visitor(), size, mapNewCoreId);

        //// updated original skyline core
        for (auto vertex: updatedVertices) {
            // find the changed skyline cores of vertex
            std::sort(vertex.newCores.begin(), vertex.newCores.end());
            std::sort(vertex.oldCores.begin(), vertex.oldCores.end());
            vector<pair<int, int>> coresAdded;
            vector<pair<int, int>> coresDeleted;
            std::set_difference(vertex.newCores.begin(), vertex.newCores.end(),
                                vertex.oldCores.begin(), vertex.oldCores.end(),
                                std::back_inserter(coresAdded));
            std::set_difference(vertex.oldCores.begin(), vertex.oldCores.end(),
                                vertex.newCores.begin(), vertex.newCores.end(),
                                std::back_inserter(coresDeleted));
            // Delete vertices in skyline cores that no longer exist
            for (int i = 0; i < coresDeleted.size(); i++) {
                int coreId = updatedMapCoCore[coresDeleted[i]];
                // remove the vertex in the old core
                auto oldEnd = std::remove(updatedSkylineCores[coreId].vertices.begin(),
                                          updatedSkylineCores[coreId].vertices.end(),
                                          vertex.vertexID);
                updatedSkylineCores[coreId].vertices.erase(oldEnd, updatedSkylineCores[coreId].vertices.end());
                if (updatedSkylineCores[coreId].vertices.size() == 0) {
                    BoundingBox bound = bounds(coresDeleted[i].first, coresDeleted[i].second, 0, 0);
                    originTree.RemoveBoundedArea(bound);
                    updatedMapCoCore.erase({coresDeleted[i].first, coresDeleted[i].second});
                }
                for (int j = 0; j < updatedSkylineCores[coreId].subGraphs.size(); j++) {
                    auto iter = std::find(updatedSkylineCores[coreId].subGraphs[j].begin(),
                                          updatedSkylineCores[coreId].subGraphs[j].end(), vertex.vertexID);
                    if (iter != updatedSkylineCores[coreId].subGraphs[j].end()) {
                        updatedSkylineCores[coreId].subGraphs[j].erase(iter);
                        if (updatedSkylineCores[coreId].subGraphs[j].size() == 0) {
                            updatedSkylineCores[coreId].subGraphs.erase(
                                    updatedSkylineCores[coreId].subGraphs.begin() + j);
                        }
                    }
                }
            }

            // Add the new vertices that appear in the skyline core
            for (int i = 0; i < coresAdded.size(); i++) {
                auto core = coresAdded[i];
                if (updatedMapCoCore.find(core) == updatedMapCoCore.end()) {
                    biSkylineCore tempSkylineCore = biSkylineCore();
                    tempSkylineCore.a = core.first;
                    tempSkylineCore.b = core.second;
                    updatedSkylineCores.push_back(tempSkylineCore);
                    updatedMapCoCore[core] = updatedSkylineCores.size() - 1;
                    originTree.Insert(updatedSkylineCores.size() - 1,
                                      bounds(tempSkylineCore.a, tempSkylineCore.b, 0, 0));
                }
                int coreId = updatedMapCoCore[core];
                updatedSkylineCores[coreId].vertices.push_back(vertex.vertexID);
                updatedSkylineCores[coreId].addedVertices.push_back(vertex.vertexID);
            }
        }
        end_time = clock();
        this_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;

        map<int, vector<pair<int, int>>> mapCoreNode;  // coreID -> nodes
        updatedResortCores(originTree, updatedMapCoNode, updatedNodeMatrix, updatedMapCoCore, updatedCoreMatrix,
                           updatedSkylineCores, mapCoreNode);

        start_time = clock();
        set<int> changedCoreId;
        set<pair<int, int>> changedNodeCo;
        //// Update the connectivity of the vertices in the core and Rtree-node
        for (int tempA = alphaMax; tempA > 0; tempA--) {
            if (tempA < alphaMax && alphaBetaMax[tempA] == alphaBetaMax[tempA + 1]) {
                bool isBlank = true;
                for (int i = 1; i <= alphaBetaMax[tempA]; i++) {
                    if (updatedCoreMatrix[tempA][i] != -1 || updatedNodeMatrix[tempA][i] != -1) {
                        isBlank = false;
                        break;
                    }
                }
                if (isBlank) {
                    continue;
                }
            }
            UnionFind UF(uNum + vNum);
            vector<bool> nodes(uNum + vNum);
            std::fill_n(nodes.begin(), nodes.size(), false);
            BoundingBox bound;
            // check if there is new deleted core
            bool start = false;
            // skip the unchanged core
//            set<int> oldVertices;
//            set<int> newVertices;
            vector<bool> newFlag(uNum + vNum);
            set<int> deletedVertices;
            for (int tempB = alphaBetaMax[tempA]; tempB > 0; tempB--) {
                vector<int> oldCores;
                vector<int> newCores;
                BoundingBox tempBound = bounds(tempA, tempB, alphaMax, 0);
//                start = clock();
                oldTree.QueryCore(RTree::AcceptEnclosing(tempBound), Visitor(), oldCores);
                newTree.QueryCore(RTree::AcceptEnclosing(tempBound), Visitor(), newCores);
//                end_time = clock();
//                total_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
                if (oldCores.size() == 0)
                    continue;
                for (auto core: newCores) {
                    for (auto vertex: newSkylineCore[core].vertices) {
                        if (!newFlag[vertex]) {
                            newFlag[vertex] = true;
                            deletedVertices.erase(vertex);
                        }
                    }
                }
//                vector<int> tempDeletedVertices;
//                for (auto it : )
                for (auto core: oldCores) {
                    for (auto vertex: oldSkylineCore[core].vertices) {
                        if (!newFlag[vertex]) {
                            deletedVertices.insert(vertex);
                        }
                    }
                }
//                if (deletedVertices.size() != 0)
//                    cout << "!!" << endl;
                if (deletedVertices.size() != 0) {
                    vector<int> tempResults;
                    //// reconstruct the subgraphs
                    vector<int> cores;
                    if (!start) {
                        bound = bounds(tempA, tempB, alphaMax, betaMax);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        //// hash table
                        for (int tempA2 = tempA; tempA2 <= alphaMax; tempA2++) {
                            if (tempB > alphaBetaMax[tempA2]) {
                                break;
                            }
                            int coreId = updatedCoreMatrix[tempA2][tempB];
                            if (coreId != -1) {
                                cores.push_back(coreId);
                            }
                        }
                    }
                    vector<pair<int, int>> tempEdges;
                    unordered_map<int, int> newVertexCore;
                    for (int j = 0; j < cores.size(); j++) {
                        int tempCoreId = cores[j];
                        vector<int> &vertices = updatedSkylineCores[tempCoreId].vertices;
                        if (vertices.size() > 0) {
                            auto it = vertices.begin();
                            int firstVertex = *it;
                            if (!nodes[firstVertex]) {
                                nodes[firstVertex] = true;
                                tempResults.push_back(firstVertex);
                            }
                            it++;
                            while (it != vertices.end()) {
//                                UF.Unite(firstVertex, *it);
                                if (!nodes[*it]) {
                                    nodes[*it] = true;
                                    tempResults.push_back(*it);
                                    newVertexCore[*it] = tempCoreId;
                                }
                                ++it;
                            }
                        }
                        for (auto edge: updatedSkylineCores[tempCoreId].connectEdges) {
                            tempEdges.push_back(edge);
                        }
                    }

//                    start_time = clock();
                    //// Calculate connectivity
                    for (int j = 0; j < tempResults.size(); j++) {
                        if (tempResults[j] < uNum) {
                            for (auto node : updatedV1Graph[tempResults[j]]) {
                                node += uNum;
                                if (nodes[node]) {
                                    UF.Unite(tempResults[j], node);
                                }
                            }
                        } else {
                            for (auto node : updatedV2Graph[tempResults[j] - uNum]) {
                                if (nodes[node]) {
                                    UF.Unite(tempResults[j], node);
                                }
                            }
                        }
                    }

                    bool change = false;
                    if (updatedCoreMatrix[tempA][tempB] != -1) {
                        int coreId = updatedCoreMatrix[tempA][tempB];
                        // Determine whether the deleted point can change the connectivity
                        for (int j = 0; j < updatedSkylineCores[coreId].subGraphs.size(); j++) {
                            map<int, vector<int>> newSubGraphs;
                            for (int &vertex : updatedSkylineCores[coreId].subGraphs[j]) {
                                newSubGraphs[UF.Find(vertex)].push_back(vertex);
                            }
                            if (newSubGraphs.size() != 1) {
                                auto position = updatedSkylineCores[coreId].subGraphs.begin() + j;
                                updatedSkylineCores[coreId].subGraphs.erase(position);
                                for (auto subgraph: newSubGraphs) {
                                    updatedSkylineCores[coreId].subGraphs.push_back(subgraph.second);
                                }
                                change = true;
                            }
                        }
                        if (change) {
                            for (auto nodeCore: mapCoreNode[coreId]) {
                                changedNodeCo.insert(nodeCore);
                            }
                        }
                    }
                    vector<pair<int, int>> newEdged;
                    for (auto &edge: tempEdges) {
                        if (deletedVertices.find(edge.first) == deletedVertices.end() &&
                            deletedVertices.find(edge.second) == deletedVertices.end()) {
                            newEdged.push_back(edge);
                        }
                    }
                } else {
                    vector<int> tempResults;
                    //// reconstruct the subgraphs
                    vector<int> cores;
                    if (!start) {
                        bound = bounds(tempA, tempB, alphaMax, betaMax);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        //// hash table
                        for (int tempA2 = tempA; tempA2 <= alphaMax; tempA2++) {
                            if (tempB > alphaBetaMax[tempA2]) {
                                break;
                            }
                            int coreId = updatedCoreMatrix[tempA2][tempB];
                            if (coreId != -1) {
                                cores.push_back(coreId);
                            }
                        }
                    }
                    vector<pair<int, int>> tempEdges;
                    unordered_map<int, int> newVertexCore;
                    for (int j = 0; j < cores.size(); j++) {
                        int tempCoreId = cores[j];
                        vector<int> &vertices = updatedSkylineCores[tempCoreId].vertices;
                        if (vertices.size() > 0) {
                            auto it = vertices.begin();
                            int firstVertex = *it;
                            if (!nodes[firstVertex]) {
                                nodes[firstVertex] = true;
                                tempResults.push_back(firstVertex);
                            }
                            it++;
                            while (it != vertices.end()) {
//                                UF.Unite(firstVertex, *it);
                                if (!nodes[*it]) {
                                    nodes[*it] = true;
                                    tempResults.push_back(*it);
                                }
                                ++it;
                            }
                        }
                    }

//                    start_time = clock();
                    //// Calculate connectivity
                    for (int j = 0; j < tempResults.size(); j++) {
                        if (tempResults[j] < uNum) {
                            for (auto node : updatedV1Graph[tempResults[j]]) {
                                node += uNum;
                                if (nodes[node]) {
                                    UF.Unite(tempResults[j], node);
                                }
                            }
                        } else {
                            for (auto node : updatedV2Graph[tempResults[j] - uNum]) {
                                if (nodes[node]) {
                                    UF.Unite(tempResults[j], node);
                                }
                            }
                        }
                    }
//                    end_time = clock();
//                    total_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
//                    if (newVertices.size() == tempCase.updatedVertices.size()) {
//                        break;
//                    }
                }

//                 updated R-tree node
                if (changedNodeCo.find({tempA, tempB}) != changedNodeCo.end()) {
                    auto RNodes = updatedMapCoNode.find({tempA, tempB})->second;
                    vector<bool> tempNodes;
                    tempNodes.resize(uNum + vNum);
                    vector<int> subgraphRoot;
                    vector<int> tempTag;
                    for (auto RNode: RNodes) {
                        int coreLeftId = RNode->left;
                        int coreRightId = RNode->right;
                        // get vertices in the range of R-Node and get the tag
                        for (int j = coreLeftId; j <= coreRightId; j++) {
                            if (j >= updatedSkylineCores.size())
                                break;
                            for (auto &tempSubgraph: updatedSkylineCores[j].subGraphs) {
                                for (auto &tempVertex: tempSubgraph) {
                                    tempNodes[tempVertex] = true;
                                }
                                int tag = -1;
                                int tempRoot = UF.Find(tempSubgraph[0]);
                                for (int k = 0; k < subgraphRoot.size(); k++) {
                                    if (subgraphRoot[k] == tempRoot) {
                                        tag = k;
                                        break;
                                    }
                                }
                                if (tag == -1) {
                                    subgraphRoot.push_back(tempRoot);
                                    tag = subgraphRoot.size() - 1;
                                }
                                tempTag.push_back(tag);
                            }
                        }
                        RNode->tag = tempTag;
                        RNode->outEdges.clear();
                        // get out edges of R-Node
                        for (int j = coreLeftId; j <= coreRightId; j++) {
                            if (j >= updatedSkylineCores.size())
                                break;
                            for (auto &tempEdge: updatedSkylineCores[j].connectEdges) {
                                if (!tempNodes[tempEdge.first] || !tempNodes[tempEdge.second]) {
                                    RNode->outEdges.push_back(tempEdge);
                                }
                            }
                        }
                    }
                    changedNodeCo.erase({tempA, tempB});
                }
            }
        }
        end_time = clock();
        this_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
        cout << this_time << "s" << endl;
        total_time += this_time;
    }
    cout << "Update time : " << total_time / (double) (cnt-1) << "s" << endl;
//    cout << "Update time : " << total_time / (double) valid << endl;
}

void deleteOneMatrix(RTree tree, BiGraph &g, int deleteNum) {
    //// load the file
    string deleteOnePath = "./ab-core/" + dataset + "/" + dataset + "-delete";
    deleteOnePath += to_string(deleteNum);
    ifstream deleteOneFile(deleteOnePath);
    string line;
    vector<deleteOneEdge> testCase;
    while (std::getline(deleteOneFile, line)) {
        if (line.size() <= 0)
            break;
        if (line[0] == 'D') {
            int pos = line.find(":");
            line = line.substr(pos + 2);
            istringstream iss(line);
            int fromVertex, toVertex;
            deleteOneEdge tempDelete;
            while (iss >> fromVertex >> toVertex) {
//                cout << fromVertex << " " << toVertex << endl;
                tempDelete.vertexPairs.push_back({fromVertex, toVertex});
            }
//            cout << tempDelete.vertexPairs.size() << endl;
            testCase.push_back(tempDelete);
        } else {
            string delimiter1 = ":";
            string delimiter2 = "|";
            int pos1 = line.find(delimiter1);
            int pos2 = line.find(delimiter2);
            string strID = line.substr(0, pos1);
            string strCoordinatesNew = line.substr(pos1 + 1, pos2 - pos1 - 1);
            string strCoordinatesOld = line.substr(pos2 + 1);

            UpdatedVertex tempVertex;
            tempVertex.vertexID = stoi(strID);

            string delimiter = ")(";

            //// Process the new coordinates
            if (strCoordinatesNew.size() > 0) {
                strCoordinatesNew = strCoordinatesNew.substr(1, strCoordinatesNew.length() - 2);
                size_t pos = 0;
                std::string token;
                while ((pos = strCoordinatesNew.find(delimiter)) != std::string::npos) {
                    token = strCoordinatesNew.substr(0, pos); // Extract a coordinate pair string
                    // Extract the x and y strings
                    size_t commaPos = token.find(",");
                    std::string x_str = token.substr(0, commaPos);
                    std::string y_str = token.substr(commaPos + 1);
                    // Convert strings to integers
                    int x = std::stoi(x_str);
                    int y = std::stoi(y_str);
                    tempVertex.newCores.push_back({x, y});
                    strCoordinatesNew.erase(0, pos + delimiter.length()); // Remove the processed part
                }
                int commaPos = strCoordinatesNew.find(",");
                std::string x_str = strCoordinatesNew.substr(0, commaPos);
                std::string y_str = strCoordinatesNew.substr(commaPos + 1);
                int x = std::stoi(x_str);
                int y = std::stoi(y_str);
                tempVertex.newCores.push_back({x, y});
            }

            //// Process the old coordinates
            if (strCoordinatesOld.size() > 0) {
                strCoordinatesOld = strCoordinatesOld.substr(1, strCoordinatesOld.length() - 2);
                size_t pos = 0;
                std::string token;
                string x_str;
                string y_str;
                int x;
                int y;
                size_t commaPos;
                while ((pos = strCoordinatesOld.find(delimiter)) != std::string::npos) {
                    token = strCoordinatesOld.substr(0, pos); // Extract a coordinate pair string
                    // Extract the x and y strings
                    commaPos = token.find(",");
                    x_str = token.substr(0, commaPos);
                    y_str = token.substr(commaPos + 1);
                    // Convert strings to integers
                    x = std::stoi(x_str);
                    y = std::stoi(y_str);
                    tempVertex.oldCores.push_back({x, y});
                    strCoordinatesOld.erase(0, pos + delimiter.length()); // Remove the processed part
                }
                commaPos = strCoordinatesOld.find(",");
                x_str = strCoordinatesOld.substr(0, commaPos);
                y_str = strCoordinatesOld.substr(commaPos + 1);
                x = std::stoi(x_str);
                y = std::stoi(y_str);
                tempVertex.oldCores.push_back({x, y});
                testCase[testCase.size() - 1].updatedVertices.push_back(tempVertex);
            }
        }
    }
    deleteOneFile.close();

    int cnt = 1;
    int valid = 0;
    clock_t start_time, end_time;
    double total_time = 0;
    for (auto tempCase: testCase) {
        double this_time = 0;
//        cout << cnt << "/" << testCase.size() << endl;
        cnt++;
        // Get skyline cores in updated vertices
        auto &updatedVertices = tempCase.updatedVertices;
        if (updatedVertices.size() == 0)
            continue;
        valid++;
        //// build temp index for test
        auto updatedV1Graph = g.neighbor_v1;
        auto updatedV2Graph = g.neighbor_v2;
        for (auto &vertexPair: tempCase.vertexPairs) {
            updatedV1Graph[vertexPair.first].erase(
                    std::remove(updatedV1Graph[vertexPair.first].begin(),
                                updatedV1Graph[vertexPair.first].end(),
                                vertexPair.second), updatedV1Graph[vertexPair.first].end());
            updatedV2Graph[vertexPair.second].erase(
                    std::remove(updatedV2Graph[vertexPair.second].begin(), updatedV2Graph[vertexPair.second].end(),
                                vertexPair.first), updatedV2Graph[vertexPair.second].end());
        }

        auto updatedSkylineCores = skylineCores;
        map<pair<int, int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>> updatedMapCoNode;
        map<pair<int, int>, int> updatedMapCoCore = mapCoCore;
        vector<vector<int>> updatedCoreMatrix = coreMatrix;
        vector<vector<int>> updatedNodeMatrix = nodeMatrix;
        RTree originTree;
        for (int i = 0; i < updatedSkylineCores.size(); i++) {
            int xx = updatedSkylineCores[i].a;
            int yy = updatedSkylineCores[i].b;
            originTree.Insert(i, bounds(xx, yy, 0, 0));
        }

        map<pair<int, int>, vector<int>> mapOldSkylineCore;
        map<pair<int, int>, vector<int>> mapNewSkylineCore;
        for (auto &tempVertex: updatedVertices) {
            for (int i = 0; i < tempVertex.oldCores.size(); i++) {
                mapOldSkylineCore[{tempVertex.oldCores[i].first, tempVertex.oldCores[i].second}].push_back(
                        tempVertex.vertexID);
            }
            for (int i = 0; i < tempVertex.newCores.size(); i++) {
                mapNewSkylineCore[{tempVertex.newCores[i].first, tempVertex.newCores[i].second}].push_back(
                        tempVertex.vertexID);
            }
        }
        // build 2 R-tree for updated vertices
        vector<vector<int>> oldMatrix;
        vector<vector<int>> newMatrix;
        int updatedAlphaMax = 0;
        vector<int> updatedAlphaVertexNum(alphaMax + 1);
        oldMatrix.resize(alphaMax + 1);
        newMatrix.resize(alphaMax + 1);
        for (int i = 0; i < alphaMax + 1; i++) {
            oldMatrix[i].resize(alphaBetaMax[i] + 1);
            newMatrix[i].resize(alphaBetaMax[i] + 1);
            fill_n(oldMatrix[i].begin(), oldMatrix[i].size(), -1);
            fill_n(newMatrix[i].begin(), newMatrix[i].size(), -1);
        }
        vector<biSkylineCore> oldSkylineCore;
        vector<biSkylineCore> newSkylineCore;
        map<pair<int, int>, int> mapOldCoreId;
        map<pair<int, int>, int> mapNewCoreId;
        int index1 = 0;
        int index2 = 0;
        for (auto &core: mapOldSkylineCore) {
            int xx = core.first.first;  // k
            int yy = core.first.second; // l
            biSkylineCore tempSkylineCore = biSkylineCore();
            tempSkylineCore.a = core.first.first;
            tempSkylineCore.b = core.first.second;
            tempSkylineCore.vertices = core.second;
            std::sort(tempSkylineCore.vertices.begin(), tempSkylineCore.vertices.end());
            oldSkylineCore.push_back(tempSkylineCore);
            mapOldCoreId[{tempSkylineCore.a, tempSkylineCore.b}] = index1;
            oldMatrix[xx][yy] = index1;
            index1++;
        }
        for (auto &core: mapNewSkylineCore) {
            int xx = core.first.first;  // k
            int yy = core.first.second; // l
            biSkylineCore tempSkylineCore = biSkylineCore();
            tempSkylineCore.a = core.first.first;
            tempSkylineCore.b = core.first.second;
            tempSkylineCore.vertices = core.second;
            std::sort(tempSkylineCore.vertices.begin(), tempSkylineCore.vertices.end());
            newSkylineCore.push_back(tempSkylineCore);
            mapNewCoreId[{tempSkylineCore.a, tempSkylineCore.b}] = index2;
            if (oldMatrix[xx][yy] != -1) {
                if (oldSkylineCore[oldMatrix[xx][yy]].vertices != tempSkylineCore.vertices) {
                    updatedAlphaMax = max(updatedAlphaMax, xx);
                }
            } else {
                updatedAlphaMax = max(updatedAlphaMax, xx);
            }
            newMatrix[xx][yy] = index2;
            index2++;
        }
//        cout << "updated alpha max:" << updatedAlphaMax << endl;

        start_time = clock();
        set<int> tempUpdatedVertex;
        for (int tempA = alphaMax; tempA > 0; tempA--) {
            for (int tempB = alphaBetaMax[tempA]; tempB > 0; tempB--) {
                if (newMatrix[tempA][tempB] != -1) {
                    for (auto &vertex : newSkylineCore[newMatrix[tempA][tempB]].vertices) {
                        tempUpdatedVertex.insert(vertex);
                    }
                }
            }
            updatedAlphaVertexNum[tempA] = tempUpdatedVertex.size();
        }

        //// updated original skyline core
        for (auto vertex: updatedVertices) {
            // find the changed skyline cores of vertex
            std::sort(vertex.newCores.begin(), vertex.newCores.end());
            std::sort(vertex.oldCores.begin(), vertex.oldCores.end());
            vector<pair<int, int>> coresAdded;
            vector<pair<int, int>> coresDeleted;
            std::set_difference(vertex.newCores.begin(), vertex.newCores.end(),
                                vertex.oldCores.begin(), vertex.oldCores.end(),
                                std::back_inserter(coresAdded));
            std::set_difference(vertex.oldCores.begin(), vertex.oldCores.end(),
                                vertex.newCores.begin(), vertex.newCores.end(),
                                std::back_inserter(coresDeleted));
            // Delete vertices in skyline cores that no longer exist
            for (int i = 0; i < coresDeleted.size(); i++) {
                int coreId = updatedMapCoCore[coresDeleted[i]];
                // remove the vertex in the old core
                auto oldEnd = std::remove(updatedSkylineCores[coreId].vertices.begin(),
                                          updatedSkylineCores[coreId].vertices.end(),
                                          vertex.vertexID);
                updatedSkylineCores[coreId].vertices.erase(oldEnd, updatedSkylineCores[coreId].vertices.end());
                if (updatedSkylineCores[coreId].vertices.size() == 0) {
                    BoundingBox bound = bounds(coresDeleted[i].first, coresDeleted[i].second, 0, 0);
                    originTree.RemoveBoundedArea(bound);
                    updatedMapCoCore.erase({coresDeleted[i].first, coresDeleted[i].second});
                }
                for (int j = 0; j < updatedSkylineCores[coreId].subGraphs.size(); j++) {
                    auto iter = std::find(updatedSkylineCores[coreId].subGraphs[j].begin(),
                                          updatedSkylineCores[coreId].subGraphs[j].end(), vertex.vertexID);
                    if (iter != updatedSkylineCores[coreId].subGraphs[j].end()) {
                        updatedSkylineCores[coreId].subGraphs[j].erase(iter);
                        if (updatedSkylineCores[coreId].subGraphs[j].size() == 0) {
                            updatedSkylineCores[coreId].subGraphs.erase(
                                    updatedSkylineCores[coreId].subGraphs.begin() + j);
                        }
                    }
                }
            }

            // Add the new vertices that appear in the skyline core
            for (int i = 0; i < coresAdded.size(); i++) {
                auto core = coresAdded[i];
                if (updatedMapCoCore.find(core) == updatedMapCoCore.end()) {
                    biSkylineCore tempSkylineCore = biSkylineCore();
                    tempSkylineCore.a = core.first;
                    tempSkylineCore.b = core.second;
                    updatedSkylineCores.push_back(tempSkylineCore);
                    updatedMapCoCore[core] = updatedSkylineCores.size() - 1;
                    originTree.Insert(updatedSkylineCores.size() - 1,
                                      bounds(tempSkylineCore.a, tempSkylineCore.b, 0, 0));
                }
                int coreId = updatedMapCoCore[core];
                updatedSkylineCores[coreId].vertices.push_back(vertex.vertexID);
                updatedSkylineCores[coreId].addedVertices.push_back(vertex.vertexID);
            }
        }
        end_time = clock();
        this_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;

        map<int, vector<pair<int, int>>> mapCoreNode;  // coreID -> nodes
        updatedResortCores(originTree, updatedMapCoNode, updatedNodeMatrix, updatedMapCoCore, updatedCoreMatrix,
                           updatedSkylineCores, mapCoreNode);

        start_time = clock();
        set<int> changedCoreId;
        set<pair<int, int>> changedNodeCo;

        UnionFind UF(uNum + vNum);
        vector<bool> nodes(uNum + vNum);
        std::fill_n(nodes.begin(), nodes.size(), false);
        vector<int> preResult;
        vector<bool> newFlag(uNum + vNum);
        vector<int> preNew;
        //// Update the connectivity of the vertices in the core and Rtree-node
        for (int tempA = updatedAlphaMax; tempA > 0; tempA--) {
            if (tempA < alphaMax && alphaBetaMax[tempA] == alphaBetaMax[tempA + 1]) {
                bool isBlank = true;
                for (int i = 1; i <= alphaBetaMax[tempA]; i++) {
                    if (updatedCoreMatrix[tempA][i] != -1 || updatedNodeMatrix[tempA][i] != -1) {
                        isBlank = false;
                        break;
                    }
                }
                if (isBlank) {
                    continue;
                }
            }
            clock_t start_time2 = clock();
            BoundingBox bound;
            UF.Init(preResult);
            for (auto &element: preResult) {
                nodes[element] = false;
            }
            preResult.clear();
            // check if there is new deleted core
            bool start = false;
            int newCnt = 0;
            set<int> deletedVertices;
            for (auto &element : preNew) {
                newFlag[element] = false;
            }
            preNew.clear();
            clock_t end_time2 = clock();
            for (int tempB = alphaBetaMax[tempA]; tempB > 0; tempB--) {
                vector<int> oldCores;
                vector<int> newCores;
                for (int tempA2 = tempA; tempA2 <= alphaMax; tempA2++) {
                    if (tempB > alphaBetaMax[tempA2]) {
                        break;
                    }
                    int coreId = oldMatrix[tempA2][tempB];
                    if (coreId != -1) {
                        oldCores.push_back(coreId);
                    }
                    coreId = newMatrix[tempA2][tempB];
                    if (coreId != -1) {
                        newCores.push_back(coreId);
                    }
                }
//                end_time = clock();
//                total_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
                if (oldCores.size() == 0)
                    continue;
                for (auto core: newCores) {
                    for (auto vertex: newSkylineCore[core].vertices) {
                        if (!newFlag[vertex]) {
                            newFlag[vertex] = true;
                            preNew.push_back(vertex);
                            deletedVertices.erase(vertex);
                            newCnt++;
                        }
                    }
                }

                for (auto core: oldCores) {
                    for (auto vertex: oldSkylineCore[core].vertices) {
                        if (!newFlag[vertex]) {
                            deletedVertices.insert(vertex);
                        }
                    }
                }

//                std::set<int> deletedVertices;
//                std::set_difference(oldVertices.begin(), oldVertices.end(),
//                                    newVertices.begin(), newVertices.end(),
//                                    std::inserter(deletedVertices, deletedVertices.begin()));
//                oldVertices = deletedVertices;
//                newVertices.clear();
                if (deletedVertices.size() != 0) {
                    vector<int> tempResults;
                    //// reconstruct the subgraphs
                    vector<int> cores;
                    if (!start) {
                        bound = bounds(tempA, tempB, alphaMax, betaMax);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        //// hash table
                        for (int tempA2 = tempA; tempA2 <= alphaMax; tempA2++) {
                            if (tempB > alphaBetaMax[tempA2]) {
                                break;
                            }
                            int coreId = updatedCoreMatrix[tempA2][tempB];
                            if (coreId != -1) {
                                cores.push_back(coreId);
                            }
                        }
                    }
                    vector<pair<int, int>> tempEdges;
                    unordered_map<int, int> newVertexCore;
                    for (int j = 0; j < cores.size(); j++) {
                        int tempCoreId = cores[j];
                        vector<int> &vertices = updatedSkylineCores[tempCoreId].vertices;
                        auto it = vertices.begin();
                        if (vertices.size() > 0) {
                            while (it != vertices.end()) {
                                if (!nodes[*it]) {
                                    nodes[*it] = true;
                                    tempResults.push_back(*it);
                                    preResult.push_back(*it);
                                    newVertexCore[*it] = tempCoreId;
                                }
                                ++it;
                            }
                        }
                        for (auto edge: updatedSkylineCores[tempCoreId].connectEdges) {
                            tempEdges.push_back(edge);
                        }
                    }

//                    start_time = clock();
                    //// Calculate connectivity
                    for (int j = 0; j < tempResults.size(); j++) {
                        if (tempResults[j] < uNum) {
                            for (auto node : updatedV1Graph[tempResults[j]]) {
                                node += uNum;
                                if (nodes[node]) {
                                    UF.Unite(tempResults[j], node);
                                }
                            }
                        } else {
                            for (auto node : updatedV2Graph[tempResults[j] - uNum]) {
                                if (nodes[node]) {
                                    UF.Unite(tempResults[j], node);
                                }
                            }
                        }
                    }

                    bool change = false;
                    if (updatedCoreMatrix[tempA][tempB] != -1) {
                        int coreId = updatedCoreMatrix[tempA][tempB];
                        // Determine whether the deleted point can change the connectivity
                        for (int j = 0; j < updatedSkylineCores[coreId].subGraphs.size(); j++) {
                            map<int, vector<int>> newSubGraphs;
                            for (int &vertex : updatedSkylineCores[coreId].subGraphs[j]) {
                                newSubGraphs[UF.Find(vertex)].push_back(vertex);
                            }
                            if (newSubGraphs.size() != 1) {
                                auto position = updatedSkylineCores[coreId].subGraphs.begin() + j;
                                updatedSkylineCores[coreId].subGraphs.erase(position);
                                for (auto subgraph: newSubGraphs) {
                                    updatedSkylineCores[coreId].subGraphs.push_back(subgraph.second);
                                }
                                change = true;
                            }
                        }
                        if (change) {
                            for (auto nodeCore: mapCoreNode[coreId]) {
                                changedNodeCo.insert(nodeCore);
                            }
                        }
                    }
                    vector<pair<int, int>> newEdged;
                    for (auto &edge: tempEdges) {
                        if (deletedVertices.find(edge.first) == deletedVertices.end() &&
                            deletedVertices.find(edge.second) == deletedVertices.end()) {
                            newEdged.push_back(edge);
                        }
                    }
                } else {
                    vector<int> tempResults;
                    //// reconstruct the subgraphs
                    vector<int> cores;
                    if (!start) {
                        bound = bounds(tempA, tempB, alphaMax, betaMax);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        //// hash table
                        for (int tempA2 = tempA; tempA2 <= alphaMax; tempA2++) {
                            if (tempB > alphaBetaMax[tempA2]) {
                                break;
                            }
                            int coreId = updatedCoreMatrix[tempA2][tempB];
                            if (coreId != -1) {
                                cores.push_back(coreId);
                            }
                        }
                    }
                    vector<pair<int, int>> tempEdges;
                    unordered_map<int, int> newVertexCore;
                    for (int j = 0; j < cores.size(); j++) {
                        int tempCoreId = cores[j];
                        vector<int> &vertices = updatedSkylineCores[tempCoreId].vertices;
                        auto it = vertices.begin();
                        int firstVertex = *it;
                        if (!nodes[firstVertex]) {
                            nodes[firstVertex] = true;
                            preResult.push_back(*it);
                            tempResults.push_back(firstVertex);
                        }
                        it++;
                        while (it != vertices.end()) {
//                                UF.Unite(firstVertex, *it);
                            if (!nodes[*it]) {
                                nodes[*it] = true;
                                tempResults.push_back(*it);
                                preResult.push_back(*it);
                            }
                            ++it;
                        }
                    }

//                    start_time = clock();
                    //// Calculate connectivity
                    for (int j = 0; j < tempResults.size(); j++) {
                        if (tempResults[j] < uNum) {
                            for (auto node : updatedV1Graph[tempResults[j]]) {
                                node += uNum;
                                if (nodes[node]) {
                                    UF.Unite(tempResults[j], node);
                                }
                            }
                        } else {
                            for (auto node : updatedV2Graph[tempResults[j] - uNum]) {
                                if (nodes[node]) {
                                    UF.Unite(tempResults[j], node);
                                }
                            }
                        }
                    }
//                    end_time = clock();
//                    total_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
                    if (newCnt >= updatedAlphaVertexNum[tempA]) {
                        break;
                    }
                }

//                 updated R-tree node
                if (changedNodeCo.find({tempA, tempB}) != changedNodeCo.end()) {
                    auto RNodes = updatedMapCoNode.find({tempA, tempB})->second;
                    vector<bool> tempNodes;
                    tempNodes.resize(uNum + vNum);
                    vector<int> subgraphRoot;
                    vector<int> tempTag;
                    for (auto RNode: RNodes) {
                        int coreLeftId = RNode->left;
                        int coreRightId = RNode->right;
                        // get vertices in the range of R-Node and get the tag
                        for (int j = coreLeftId; j <= coreRightId; j++) {
                            if (j >= updatedSkylineCores.size())
                                break;
                            for (auto &tempSubgraph: updatedSkylineCores[j].subGraphs) {
                                for (auto &tempVertex: tempSubgraph) {
                                    tempNodes[tempVertex] = true;
                                }
                                int tag = -1;
                                int tempRoot = UF.Find(tempSubgraph[0]);
                                for (int k = 0; k < subgraphRoot.size(); k++) {
                                    if (subgraphRoot[k] == tempRoot) {
                                        tag = k;
                                        break;
                                    }
                                }
                                if (tag == -1) {
                                    subgraphRoot.push_back(tempRoot);
                                    tag = subgraphRoot.size() - 1;
                                }
                                tempTag.push_back(tag);
                            }
                        }
                        RNode->tag = tempTag;
                        RNode->outEdges.clear();
                        // get out edges of R-Node
                        for (int j = coreLeftId; j <= coreRightId; j++) {
                            if (j >= updatedSkylineCores.size())
                                break;
                            for (auto &tempEdge: updatedSkylineCores[j].connectEdges) {
                                if (!tempNodes[tempEdge.first] || !tempNodes[tempEdge.second]) {
                                    RNode->outEdges.push_back(tempEdge);
                                }
                            }
                        }
                    }
                    changedNodeCo.erase({tempA, tempB});
                }
            }
        }
        end_time = clock();
        this_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
//        cout << this_time << "s" << endl;
        total_time += this_time;
    }
    cout << "Update time : " << total_time / (double) (cnt-1) << "s" << endl;
//    cout << "Update time : " << total_time / (double) valid << endl;
}

int main(int argc, char **argv) {
    srand(time(0));
    RTree tree;
    Visitor x;

    timeval start_time, end_time;
    dataset =argv[1];
//    dataset = "TREC";
//    dataset = "Flickr";
//    dataset = "WebTrackers";
//    string dataset = "test";
    // load skyline cores
    string datasetPath = "./ab-core/" + dataset + "/" + dataset;
    BiGraph g(datasetPath);
    uNum = g.num_v1;
    vNum = g.num_v2;
    ifstream inFile;
    string filePath = "./ab-core/" + dataset + "/" + dataset + "-skylineCore";
    inFile.open(filePath);
    cout << "Loading skyline cores..." << endl;
    loadBiSkylineCores(inFile, skylineCores, alphaMax, betaMax);
    inFile.close();

    cout << "alpha Max: " << alphaMax << "beta Max: " << betaMax << endl;

//    // preprocess
    preProcess(skylineCores, g);

    cout << skylineCores.size() << endl;

    filePath = "./ab-core/" + dataset + "/" + dataset + "-alphaBetaMax";
    inFile.open(filePath);
    loadAlphaBetaMax(inFile);
    inFile.close();

    int coreNum = skylineCores.size();
    long long indexSize = 0;
    for (int i = 0; i < skylineCores.size(); i++) {
        indexSize += 4 * 2;
        indexSize += skylineCores[i].vertexU.size() * 4;
        indexSize += skylineCores[i].vertexV.size() * 4;
    }
    cout << "Index Size of Skyline Cores: " << indexSize / 1024.0 / 1024.0 << "MB" << endl;

    // build R-tree
    gettimeofday(&start_time, nullptr); // record start time
    for (int i = 0; i < coreNum; i++) {
        int xx = skylineCores[i].a;
        int yy = skylineCores[i].b;
        tree.Insert(i, bounds(xx, yy, 0, 0));
    }
    gettimeofday(&end_time, nullptr); // record end time
    double buildTime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1e6; // 计算时间间隔
    std::cout << "Building Time :" << buildTime << "s" << std::endl;

    // calculate index size of R*tree
//    x = tree.Query(RTree::AcceptAny(), Visitor());
    tree.indexSize(RTree::AcceptAny(), Visitor(), indexSize);
    resortCores(tree);
//    cout << "Index Size of R-tree and Skyline Cores: " << indexSize / 1024.0 / 1024.0 << "MB" << endl;

    auto root = tree.m_root;
    int height = 0;
    while (root != NULL) {
        height++;
        if (!root->hasLeaves) {
            root = static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(root->items[0]);
        } else {
            break;
        }
    }
//    addOne(tree, g, 1);

    gettimeofday(&start_time, nullptr); // record start time
    connectivity4(g, tree);
//    connectivity1(g, tree);
    gettimeofday(&end_time, nullptr); // record end time
    double connectTime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1e6; // 计算时间间隔
    std::cout << "connect Time :" << connectTime << "s" << std::endl;
    cout << "R-tree search time : " << RtreeTime << "s" << endl;

    cout << "Index Size of R-tree and Skyline Cores: " << indexSize / 1024.0 / 1024.0 << "MB" << endl;
    int subgraphNum = 0;
    int connectEdgeNum = 0;
    for (auto &core : skylineCores) {
        subgraphNum += core.subGraphs.size();
        connectEdgeNum += core.connectEdges.size();
    }
    indexSize += height * subgraphNum * 4;
    indexSize += (outEdgeNum + connectEdgeNum) * 8;
    cout << "R-tree Height : " << height << endl;
    cout << "The number of skyline cores : " << skylineCores.size() << endl;
    cout << "The number of subgraphs : " << subgraphNum << endl;
    cout << "The number of connect edges : " << connectEdgeNum << endl;
    cout << "The number of out edges : " << outEdgeNum << endl;

    cout << "Index Size of R-tree and Skyline Cores and Tags and Edges: " << indexSize / 1024.0 / 1024.0 << "MB" << endl;

//    vector<int> nums = {1, 4, 16, 64, 256, 1024};
    cout << "ADD TEST:" << endl;
    vector<int> nums = {1, 4, 16, 64, 256, 1024};
    for (int i = 0; i < nums.size(); i++) {
        addOneMatrix(tree, g, nums[i]);
    }
    cout << "DELETE TEST:" << endl;
    for (int i = 0; i < nums.size(); i++) {
        deleteOneMatrix(tree, g, nums[i]);
    }

    std::cout << "connect Time :" << connectTime << "s" << std::endl;

    return 0;
}



