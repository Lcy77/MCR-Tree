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

ofstream expData;
int pointNum;

BoundingBox bounds(int x, int y, int w, int h) {
    BoundingBox bb;

    bb.edges[0].first = x;
    bb.edges[0].second = x + w;

    bb.edges[1].first = y;
    bb.edges[1].second = y + h;

    return bb;
}

struct Visitor {
    int count;
    bool ContinueVisiting;

    Visitor() : count(0), ContinueVisiting(true) {};

    void operator()(const RTree::Leaf *const leaf) {
        count++;
    }
};

void getSubGraphUF(vector<int> &cores, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *> &nodes,
                   BiGraph &g,  vector<bool> &isResult, UnionFind &UF) {
    clock_t start_time, end_time;
    subBiGraph tempSubgraph;
    vector<int> results;
    vector<pair<int, int>> edges;

    //// process cores
    for (int i = 0; i < cores.size(); i++) {
        int coreId = cores[i];
        for (int j = 0; j < skylineCores[coreId].subGraphs.size(); j++) {
            vector<int> &vertices = skylineCores[coreId].subGraphs[j];
            if (vertices.size() > 0) {
                auto it = vertices.begin();
                int firstVertex = *it;
                if (!isResult[firstVertex]) {
                    isResult[firstVertex] = true;
                    results.push_back(firstVertex);
                }
                it++;
                while (it != vertices.end()) {
                    UF.Unite(firstVertex, *it);
                    if (!isResult[*it]) {
                        isResult[*it] = true;
                        results.push_back(*it);
                    }
                    ++it;
                }
            }
        }
        for (auto &tempEdge: skylineCores[coreId].connectEdges) {
            edges.push_back(tempEdge);
        }
    }

    //// process r-tree nodes
    for (auto &node: nodes) {
        int left = node->left;
        int right = node->right;
        unordered_map<int, int> connectVertex;
        for (int i = left; i <= right; i++) {
            int cnt = 0;
            for (int j = 0; j < skylineCores[i].subGraphs.size(); j++) {
                vector<int> &vertices = skylineCores[i].subGraphs[j];
                auto it = vertices.begin();
                int firstVertex = *it;
                if (!isResult[firstVertex]) {
                    isResult[firstVertex] = true;
                    results.push_back(firstVertex);
                }
                it++;
                while (it != vertices.end()) {
                    UF.Unite(firstVertex, *it);
                    if (!isResult[*it]) {
                        isResult[*it] = true;
                        results.push_back(*it);
                    }
                    ++it;
                }
                if (connectVertex.find(cnt) == connectVertex.end()) {
                    connectVertex[cnt] = firstVertex;
                } else {
                    UF.Unite(connectVertex[cnt], vertices[0]);
                }
                cnt++;
            }
        }
        for (auto &tempEdge : node->outEdges) {
            edges.push_back(tempEdge);
        }
    }
    for (auto &tempEdge : edges) {
        if (isResult[tempEdge.first] && isResult[tempEdge.second]) {
            UF.Unite(tempEdge.first, tempEdge.second);
        }
    }
//    for (auto vertex : results) {
//        if (vertex < g.num_v1) {
//            for (auto neighbor : g.neighbor_v1[vertex]) {
//                if (isResult[neighbor + g.num_v1]) {
//                    UF.Unite(vertex, neighbor + g.num_v1);
//                }
//            }
//        } else {
//            for (auto neighbor : g.neighbor_v2[vertex - g.num_v1]) {
//                if (isResult[neighbor]) {
//                    UF.Unite(vertex, neighbor);
//                }
//            }
//        }
//    }
    map<int, vector<int>> rootGraph;
    for (auto &item : results) {
        rootGraph[UF.Find(item)].push_back(item);
    }
    cout << rootGraph.size() << endl;
//    cout << results.size() << " ";
}


void queryTest(RTree &tree, BiGraph &g) {
//    int queryAB = 1;
//    if (dataset == "TREC")
//        queryAB = 509;
//    else if (dataset == "Flickr")
//        queryAB = 148;
    vector<bool> left_node;
    vector<bool> right_node;
    int testCase = 4;
    int cnt = 0;
//    timeval start_time, end_time;
    clock_t start_time, end_time;
    double totalQueryTime;

    vector<pair<int, int>> list = {{123, 80}, {70, 187}, {122, 183}, {106, 180}};
    for (int i = 0; i < testCase; i++) {
        printf("%d / %d   ", i + 1, testCase);
        alpha = rand() % alphaMax;
//      beta = rand() % betaMax / 10 + 1;
        beta = rand() % alphaBetaMax[alpha] + 1;
        alpha = list[i].first;
        beta = list[i].second;
//        alpha = rand() % 1000 + 1;
//        beta = rand() % 1000 + 1;
//        alpha = queryAB / 2.0;
//        beta = queryAB / 2.0;
        cout << "alpha:" << alpha << " " << "beta:" << beta << "  ";
        BoundingBox bound = bounds(alpha, beta, alphaMax, betaMax);
        vector<int> cores;
        vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *> nodes;
        double queryTime = 0.0;
        start_time = clock();
        tree.QueryPrune(RTree::AcceptEnclosing(bound), Visitor(), cores, nodes);
        end_time = clock();
        queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
        if (cores.size() > 0 || nodes.size() > 0) {
            vector<bool> isResult;
            isResult.resize(uNum + vNum);
            std::fill_n(isResult.begin(), isResult.size(), false);
            UnionFind UF(uNum + vNum);
            start_time = clock();
            getSubGraphUF(cores, nodes, g, isResult, UF);
            end_time = clock();
            queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;

            cout << queryTime << endl;
            totalQueryTime += queryTime;
            cnt++;
        } else {
            cout << endl;
        }
    }
    cout << "Query time: " << totalQueryTime / (double) cnt << "ms" << endl;
}

void queryTestFixAlpha(RTree &tree, BiGraph &g, int fixedAlpha, vector<int> betaNums) {
    ofstream resultFile;
    string resultFilePath = "./ab-core/" + dataset + "/" + dataset + "_query_a";
    resultFile.open(resultFilePath);
    vector<bool> left_node;
    vector<bool> right_node;
    int cnt = 0;
//    timeval start_time, end_time;
    clock_t start_time, end_time;

    for (int i = 0; i < betaNums.size(); i++) {
        alpha = fixedAlpha;
        beta = betaNums[i];
        cout << "alpha:" << alpha << " " << "beta:" << beta << "  ";
        resultFile << "alpha:" << alpha << " " << "beta:" << beta << "  ";
        BoundingBox bound = bounds(alpha, beta, alphaMax, betaMax);
        vector<int> cores;
        vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *> nodes;
        double queryTime = 0.0;
        start_time = clock();
        tree.QueryPrune(RTree::AcceptEnclosing(bound), Visitor(), cores, nodes);
        end_time = clock();
        queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
        if (cores.size() > 0) {
            vector<bool> isResult;
            isResult.resize(uNum + vNum);
            std::fill_n(isResult.begin(), isResult.size(), false);
            UnionFind UF(uNum + vNum);
            start_time = clock();
            getSubGraphUF(cores, nodes, g, isResult, UF);
            end_time = clock();
            queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
            cout << queryTime << endl;
            resultFile << queryTime << endl;
            cnt++;
        } else {
            cout << endl;
        }
    }
}

void queryTestFixBeta(RTree &tree, BiGraph &g, int fixedBeta, vector<int> alphaNums) {
    ofstream resultFile;
    string resultFilePath = "./ab-core/" + dataset + "/" + dataset + "_query_b";
    resultFile.open(resultFilePath);
    vector<bool> left_node;
    vector<bool> right_node;
    int cnt = 0;
//    timeval start_time, end_time;
    clock_t start_time, end_time;

    for (int i = 0; i < alphaNums.size(); i++) {
        beta = fixedBeta;
        alpha = alphaNums[i];
        cout << "alpha:" << alpha << " " << "beta:" << beta << "  ";
        resultFile << "alpha:" << alpha << " " << "beta:" << beta << "  ";
        BoundingBox bound = bounds(alpha, beta, alphaMax, betaMax);
        vector<int> cores;
        vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *> nodes;
        double queryTime = 0.0;
        start_time = clock();
        tree.QueryPrune(RTree::AcceptEnclosing(bound), Visitor(), cores, nodes);
        end_time = clock();
        queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
        if (cores.size() > 0) {
            vector<bool> isResult;
            isResult.resize(uNum + vNum);
            std::fill_n(isResult.begin(), isResult.size(), false);
            UnionFind UF(uNum + vNum);
            start_time = clock();
            getSubGraphUF(cores, nodes, g, isResult, UF);
            end_time = clock();
            queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
            cout << queryTime << endl;
            resultFile << queryTime << endl;
            cnt++;
        } else {
            cout << endl;
        }
    }
    resultFile.close();
}

void queryTestRand(RTree &tree, BiGraph &g) {
    vector<bool> left_node;
    vector<bool> right_node;
    int cnt = 0;
//    timeval start_time, end_time;
    clock_t start_time, end_time;
    int testCase = 10000;
    double total_time = 0;

    for (int i = 0; i < testCase; i++) {
        alpha = rand() % (alphaMax - 1) + 1;
        if (alphaBetaMax[alpha] == 1)
            beta = 1;
        else
            beta = rand() % (alphaBetaMax[alpha] - 1) + 1;
        cout << "alpha:" << alpha << " " << "beta:" << beta << "  ";
        BoundingBox bound = bounds(alpha, beta, alphaMax, betaMax);
        vector<int> cores;
        vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *> nodes;
        double queryTime = 0.0;
        start_time = clock();
        tree.QueryPrune(RTree::AcceptEnclosing(bound), Visitor(), cores, nodes);
        end_time = clock();
        queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
        if (cores.size() > 0) {
            vector<bool> isResult;
            isResult.resize(uNum + vNum);
            std::fill_n(isResult.begin(), isResult.size(), false);
            UnionFind UF(uNum + vNum);
            start_time = clock();
            getSubGraphUF(cores, nodes, g, isResult, UF);
            end_time = clock();
            queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
            cout << queryTime << endl;
            cnt++;
        } else {
            cout << endl;
        }
        total_time += queryTime;
    }
    cout << "query time : " << total_time / (double) cnt << "ms" << endl;
    expData << "query time : " << total_time / (double) cnt << "ms" << endl;
}


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
                            int tempCoreId = newVertexCore[tempResults[j]];
                            skylineCores[tempCoreId].connectEdges.push_back({tempResults[j], node});
                        }
                    }
                } else {
                    for (int node: g.neighbor_v2[tempResults[j] - uNum]) {
                        if (nodes[node] && UF2.Find(tempResults[j]) != UF2.Find(node)) {
                            UF2.Unite(tempResults[j], node);
                            int tempCoreId = newVertexCore[tempResults[j]];
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
    string expDataPath = "./ab-core/" + dataset + "/" + dataset + "-expData.txt";
    expData.open(expDataPath);

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
    expData << "Building Time :" << buildTime << "s" << std::endl;

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

    gettimeofday(&start_time, nullptr); // record start time
    connectivity4(g, tree);
//    connectivity1(g, tree);
    gettimeofday(&end_time, nullptr); // record end time
    double connectTime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1e6; // 计算时间间隔
    std::cout << "connect Time :" << connectTime << "s" << std::endl;
    expData << "connect Time :" << connectTime << "s" << std::endl;
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

    expData << "R-tree Height : " << height << endl;
    expData << "The number of skyline cores : " << skylineCores.size() << endl;
    expData << "The number of subgraphs : " << subgraphNum << endl;
    expData << "The number of connect edges : " << connectEdgeNum << endl;
    expData << "The number of out edges : " << outEdgeNum << endl;

    expData << "Index Size of R-tree and Skyline Cores and Tags and Edges: " << indexSize / 1024.0 / 1024.0 << "MB" << endl;

//    queryTest(tree, g);

    if (dataset == "Flickr") {
        int fixedAlpha = 140;
        vector<int> betaNums = {5, 20, 35, 50, 65, 80, 95, 110, 125, 140};
        queryTestFixAlpha(tree, g, fixedAlpha, betaNums);
        int fixedBeta = 140;
        vector<int> alphaNums = {5, 20, 35, 50, 65, 80, 95, 110, 125, 140};
        queryTestFixBeta(tree, g, fixedBeta, alphaNums);
    } else if (dataset == "WebTrackers" || dataset == "TREC" || dataset == "Orkut") {
        int fixedAlpha = 400;
        vector<int> betaNums = {40, 80, 120, 160, 200, 240, 280, 320, 360, 400};
        queryTestFixAlpha(tree, g, fixedAlpha, betaNums);
        int fixedBeta = 400;
        vector<int> alphaNums = {40, 80, 120, 160, 200, 240, 280, 320, 360, 400};
        queryTestFixBeta(tree, g, fixedBeta, alphaNums);
    }

    queryTestRand(tree, g);
//    queryTestRand2(tree, g);

    expData.close();

    return 0;
}



