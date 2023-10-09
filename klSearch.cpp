#include <string>
#include <ctime>
#include <sys/time.h>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <stdio.h>
#include <queue>
#include <deque>
#include <omp.h>
#include "RStarTree.h"
#include "kl-core/DiGraph.h"
#include "kl-core/KlSkylineCore.h"
#include "utilities/UnionFind.h"

using namespace std;
const int dimensions = 2;
const int min_child_items = 32;
const int max_child_items = 64;
typedef RStarTree<int, dimensions, min_child_items, max_child_items> RTree;

typedef RTree::BoundingBox BoundingBox;

int kMax, lMax;
int k, l;
int vertexNum;
vector<KlSkylineCore> skylineCores;
map<int, vector<int>> mapVertexCore;
map<pair<int, int>, int> mapCoreId;
vector<int> klMax;
map<pair<int, int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>> mapCoNode;
map<pair<int, int>, int> mapCoCore;
//map<pair<int, int>, int> mapConnectNum;
vector<vector<int>> coreMatrix;
vector<vector<int>> nodeMatrix;

//ofstream expData;

int pointNum;

//string dataset = "BerkStan";
//    string dataset = "wiki";
//    string dataset = "arabic-2005";
//    string dataset = "in-2004";

string dataset = "";

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
              DiGraph &g,  vector<bool> &isResult, UnionFind &UF) {
    clock_t start_time, end_time;
    subBiGraph tempSubgraph;

    vector<int> results;
    vector<pair<int, int>> edges;
    // init the Union Find

    start_time = clock();
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
    end_time = clock();

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
    map<int, vector<int>> rootMap;
    for (auto vertex : results) {
        rootMap[UF.Find(vertex)].push_back(vertex);
    }
}


void queryTest(RTree &tree, DiGraph &g) {
    int testCase = 1;
    int cnt = 0;
//    timeval start_time, end_time;
    clock_t start_time, end_time;
    double totalQueryTime;

//    vector<pair<int, int>> list = {{42,43}, {11, 21}, {38, 23}};
    for (int i = 0; i < testCase; i++) {
        printf("%d / %d   ", i + 1, testCase);
        k = 0;
//      beta = rand() % betaMax / 10 + 1;
        l = 1;
        cout << "k:" << k << " " << "l:" << l << "  ";
        BoundingBox bound = bounds(k, l, kMax, lMax);
        vector<int> cores;
        vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *> nodes;
        double queryTime = 0.0;
        start_time = clock();
//        tree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
        tree.QueryPrune(RTree::AcceptEnclosing(bound), Visitor(), cores, nodes);
//        cout << "core number : " << cores.size() << endl;
//        cout << "node number : " << nodes.size() << endl;
        end_time = clock();
        queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
        if (cores.size() > 0 || nodes.size() > 0) {
            vector<bool> isResult;
            isResult.resize(vertexNum);
            std::fill_n(isResult.begin(), isResult.size(), false);
            UnionFind UF(vertexNum);
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

void queryTestFixK(RTree &tree, DiGraph &g, int fixedK, vector<int> lNums) {
    ofstream resultFile;
    string resultFilePath = "./kl-core/" + dataset + "/" + dataset + "_query_k";
    resultFile.open(resultFilePath);
//    timeval start_time, end_time;
    clock_t start_time, end_time;

//    vector<pair<int, int>> list = {{42,43}, {11, 21}, {38, 23}};
    for (int i = 0; i < lNums.size(); i++) {
        k = fixedK;
//      beta = rand() % betaMax / 10 + 1;
        l = lNums[i];
        cout << "k:" << k << " " << "l:" << l << "  ";
        resultFile << "k:" << k << " " << "l:" << l << "  ";
        BoundingBox bound = bounds(k, l, kMax, lMax);
        vector<int> cores;
        vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *> nodes;
        double queryTime = 0.0;
        start_time = clock();
//        tree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
        tree.QueryPrune(RTree::AcceptEnclosing(bound), Visitor(), cores, nodes);
//        cout << "core number : " << cores.size() << endl;
//        cout << "node number : " << nodes.size() << endl;
        end_time = clock();
        queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
        if (cores.size() > 0 || nodes.size() > 0) {
            vector<bool> isResult;
            isResult.resize(vertexNum);
            std::fill_n(isResult.begin(), isResult.size(), false);
            UnionFind UF(vertexNum);
            start_time = clock();
            getSubGraphUF(cores, nodes, g, isResult, UF);
            end_time = clock();
            queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
            cout << queryTime << endl;
            resultFile << queryTime << endl;
        } else {
            cout << endl;
        }
    }
    resultFile.close();
}

void queryTestFixL(RTree &tree, DiGraph &g, int fixedL, vector<int> kNums) {
    ofstream resultFile;
    string resultFilePath = "./kl-core/" + dataset + "/" + dataset + "_query_l";
    resultFile.open(resultFilePath);
//    timeval start_time, end_time;
    clock_t start_time, end_time;

//    vector<pair<int, int>> list = {{42,43}, {11, 21}, {38, 23}};
    for (int i = 0; i < kNums.size(); i++) {
        l = fixedL;
//      beta = rand() % betaMax / 10 + 1;
        k = kNums[i];
        cout << "k:" << k << " " << "l:" << l << "  ";
        resultFile << "k:" << k << " " << "l:" << l << "  ";
        BoundingBox bound = bounds(k, l, kMax, lMax);
        vector<int> cores;
        vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *> nodes;
        double queryTime = 0.0;
        start_time = clock();
//        tree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
        tree.QueryPrune(RTree::AcceptEnclosing(bound), Visitor(), cores, nodes);
//        cout << "core number : " << cores.size() << endl;
//        cout << "node number : " << nodes.size() << endl;
        end_time = clock();
        queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
        if (cores.size() > 0 || nodes.size() > 0) {
            vector<bool> isResult;
            isResult.resize(vertexNum);
            std::fill_n(isResult.begin(), isResult.size(), false);
            UnionFind UF(vertexNum);
            start_time = clock();
            getSubGraphUF(cores, nodes, g, isResult, UF);
            end_time = clock();
            queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
            resultFile << queryTime << endl;
            cout << queryTime << endl;
        } else {
            cout << endl;
        }
    }
}


bool isEqual(std::vector<int> const &a, std::vector<int> const &b) {
    if (a.size() != b.size()) return false;
    return std::equal(a.begin(), a.end(), b.begin());
}

void preProcess(vector<KlSkylineCore> &skylineCores, DiGraph &g) {
    int cnt = 0;
    for (int i = 0; i < skylineCores.size(); i++) {
        mapCoreId.insert(make_pair(make_pair(skylineCores[i].k, skylineCores[i].l), i));
        for (int j = 0; j < skylineCores[i].vertices.size(); j++) {
            int node = skylineCores[i].vertices[j];
            mapVertexCore[node].push_back(i);
        }
        printf("%d / %d\n", i, skylineCores.size());
    }
}

void loadKLMax(ifstream &ifstream1) {
    klMax.resize(kMax + 1);
    coreMatrix.resize(kMax + 1);
    nodeMatrix.resize(kMax + 1);
    //// load alpha beta max file
    for (int i = 0; i <= kMax; i++) {
        int tempK;
        int tempL;
        ifstream1 >> tempK >> tempL;
        klMax[tempK] = tempL;
        pointNum += tempL;
        coreMatrix[tempK].resize(tempL + 1);
        nodeMatrix[tempK].resize(tempL + 1);
        std::fill(coreMatrix[tempK].begin(), coreMatrix[tempK].end(), -1);
        std::fill(nodeMatrix[tempK].begin(), nodeMatrix[tempK].end(), -1);
    }
}

void rtreeDFS(RStarTree<int, dimensions, min_child_items, max_child_items>::Node *node, int &coreCnt,
              vector<KlSkylineCore> &newSkylineCores) {
    mapCoNode[{node->bound.edges[0].first, node->bound.edges[1].first}].push_back(node);
    nodeMatrix[node->bound.edges[0].first][node->bound.edges[1].first] = 1;
    if (node->hasLeaves) {
        node->left = coreCnt;
        for (auto item1: node->items) {
            static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Leaf *>(item1)->coreId = coreCnt;
            newSkylineCores.push_back(
                    skylineCores[mapCoreId[{item1->bound.edges[0].first, item1->bound.edges[1].first}]]);
            newSkylineCores.back().id = coreCnt;
            mapCoCore[{item1->bound.edges[0].first, item1->bound.edges[1].first}] = coreCnt;
            coreCnt++;
        }
        node->right = coreCnt - 1;
    } else {
        for (auto item1: node->items) {
            rtreeDFS(static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(item1), coreCnt,
                     newSkylineCores);
        }
        node->left = static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(node->items[0])->left;
        node->right = static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(node->items.back())->right;
//        cout << node->left << " " << node->right << endl;
    }
}

void resortCores(RTree &tree) {
    auto node = tree.m_root;
    int coreCnt = 0;
    vector<KlSkylineCore> newSkylineCores;
//    map<pair<int, int>, int> newMapCoreId;
    if (node != NULL) {
        rtreeDFS(node, coreCnt, newSkylineCores);
    }
    skylineCores = newSkylineCores;
    for (int i = 0; i < skylineCores.size(); i++) {
        coreMatrix[skylineCores[i].k][skylineCores[i].l] = i;
    }
//    for (auto &item : mapCoNode) {
//        for (auto node : item.second) {
//            cout << "range:" << node->left << " " << node->right << endl;
//        }
//    }
}

bool compareSkylineCore(const KlSkylineCore &a, const KlSkylineCore &b) {
    if (a.l != b.l) {
        return a.l > b.l;
    } else {
        return a.k < b.k;
    }
}

// add Tags and Connected Edges
int outEdgeNum = 0;

void connectivity(DiGraph &g, RTree &tree) {
    timeval start_time, end_time;
    int cnt = 0;
    UnionFind UF(vertexNum);
    UnionFind UF2(vertexNum);
    vector<bool> nodes(vertexNum);
    std::fill_n(nodes.begin(), nodes.size(), false);
    vector<int> preResult;
    for (int tempK = kMax; tempK >= 0; tempK--) {
        cout << "k:" << tempK << "/" << kMax << endl;
        if (tempK < kMax && klMax[tempK] == klMax[tempK + 1]) {
            bool isBlank = true;
            for (int i = 1; i <= klMax[tempK]; i++) {
//                cout << nodeMatrix[a][i] << endl;
                if (coreMatrix[tempK][i] != -1 || nodeMatrix[tempK][i] != -1) {
                    isBlank = false;
                    break;
                }
            }
//            cout << isBlank << endl;
            if (isBlank) {
                cnt += klMax[tempK];
                continue;
            }
        }
        UF.Init(preResult);
        UF2.Init(preResult);
        for (auto &element: preResult) {
            nodes[element] = false;
        }
        preResult.clear();
//        BoundingBox bound;
//        vector<pair<int, int>> preCriticalEdges;
        for (int tempL = klMax[tempK]; tempL >= 0; tempL--) {
            if (cnt % 10000 == 0)
                cout << cnt << "/" << pointNum << endl;
            cnt++;
            vector<int> cores;  // dominate core for the current core
            vector<int> tempResults;
//            bound = bounds(a, b, alphaMax, 0);
//            tree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
            //// hash table
            for (int tempA = tempK; tempA <= kMax; tempA++) {
                if (tempL > klMax[tempA]) {
                    break;
                }
                int coreId = coreMatrix[tempA][tempL];
                if (coreId != -1) {
                    cores.push_back(coreId);
                }
            }
            unordered_map<int, int> newVertexCore;
            //// get the newly added vertices
            for (int j = 0; j < cores.size(); j++) {
                int coreId = cores[j];
                vector<int> &vertices = skylineCores[coreId].vertices;
                auto it = vertices.begin();
                while (it != vertices.end()) {
                    if (!nodes[*it]) {
                        nodes[*it] = true;
                        preResult.push_back(*it);
                        tempResults.push_back(*it);
                        newVertexCore[*it] = coreId;
                    }
                    ++it;
                }
            }

            //// Calculate connectivity
            for (int j = 0; j < tempResults.size(); j++) {
                for (int node: g.inGraph[tempResults[j]]) {
                    if (nodes[node] && UF.Find(tempResults[j]) != UF.Find(node)) {
                        UF.Unite(tempResults[j], node);
                    }
                }
                for (int node: g.outGraph[tempResults[j]]) {
                    if (nodes[node] && UF.Find(tempResults[j]) != UF.Find(node)) {
                        UF.Unite(tempResults[j], node);
                    }
                }
            }

            //// Count connected subgraphs of current core
            if (coreMatrix[tempK][tempL] != -1) {
                int coreID = coreMatrix[tempK][tempL];
                map<int, vector<int>> rootVec;
                for (auto vertex: skylineCores[coreID].vertices) {
                    int tempRoot = UF.Find(vertex);
                    rootVec[tempRoot].push_back(vertex);
                }
                for (auto &item1: rootVec) {
                    skylineCores[coreID].subGraphs.push_back(item1.second);
                }
            }

            //// Count the connected edges
            for (int j = 0; j < cores.size(); j++) {
                int coreId = cores[j];
                for (auto &subGraph: skylineCores[coreId].subGraphs) {
                    auto it = subGraph.begin();
                    int firstVertex = *it;
                    it++;
                    while (it != subGraph.end()) {
                        UF2.Unite(firstVertex, *it);
                        it++;
                    }
                }
                for (auto &edge: skylineCores[coreId].connectEdges) {
                    UF2.Unite(edge.first, edge.second);
                }
            }
            for (int j = 0; j < tempResults.size(); j++) {
                for (int node: g.inGraph[tempResults[j]]) {
                    if (nodes[node] && UF2.Find(tempResults[j]) != UF2.Find(node)) {
                        UF2.Unite(tempResults[j], node);
                        int tempCoreId = newVertexCore[tempResults[j]];
                        skylineCores[tempCoreId].connectEdges.push_back({tempResults[j], node});
                    }
                }
                for (int node: g.outGraph[tempResults[j]]) {
                    if (nodes[node] && UF2.Find(tempResults[j]) != UF2.Find(node)) {
                        UF2.Unite(tempResults[j], node);
                        int tempCoreId = newVertexCore[tempResults[j]];
                        skylineCores[tempCoreId].connectEdges.push_back({tempResults[j], node});
                    }
                }
            }

            //// calculate tags and out edges for R-tree node
            if (mapCoNode.find({tempK, tempL}) != mapCoNode.end()) {
                auto RNodes = mapCoNode.find({tempK, tempL})->second;
                vector<bool> tempNodes;
                tempNodes.resize(vertexNum);
                vector<int> subgraphRoot;
                vector<int> tempTag;
                for (auto RNode: RNodes) {
                    int coreLeftId = RNode->left;
                    int coreRightId = RNode->right;
                    // get vertices in the range of R-Node and get the tag
                    for (int j = coreLeftId; j <= coreRightId; j++) {
                        for (auto &tempSubgraph: skylineCores[j].subGraphs) {
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
                    for (int j = coreLeftId; j <= coreRightId; j++) {
                        for (auto &tempEdge: skylineCores[j].connectEdges) {
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
    dataset = argv[1];
    string datasetPath = "./kl-core/" + dataset + "/" + dataset;

    string expDataFile = "./kl-core/" + dataset + "/" + dataset + "-expData.txt";
//    expData.open(expDataFile);

    // load skyline cores
    ifstream inFile;
    string skylinePath = "./kl-core/" + dataset + "/" + dataset + "_skyline";
    inFile.open(skylinePath);
    cout << "Loading skyline cores..." << endl;
    loadKlSkylineCores(inFile, skylineCores, kMax, lMax);
    inFile.close();

    string degreeFile = datasetPath + "_degree.dat";
    string graphFile = datasetPath + "_graph.dat";
    // load graph
    DiGraph g(degreeFile, graphFile);
    g.read_graph();
    vertexNum = g.numOfNodes;

    // preprocess
    preProcess(skylineCores, g);

    cout << skylineCores.size() << endl;

    string filePath = "./kl-core/" + dataset + "/" + dataset + "_klMax";
    inFile.open(filePath);
    loadKLMax(inFile);
    inFile.close();

    int coreNum = skylineCores.size();
    long long indexSize = 0;
    for (int i = 0; i < skylineCores.size(); i++) {
        indexSize += 4 * 2;
        indexSize += skylineCores[i].vertices.size() * 4;
    }
    cout << "Index Size of Skyline Cores: " << indexSize / 1024.0 / 1024.0 << "MB" << endl;

    // build R-tree
    gettimeofday(&start_time, nullptr); // record start time
    for (int i = 0; i < coreNum; i++) {
        int xx = skylineCores[i].k;
        int yy = skylineCores[i].l;
        tree.Insert(i, bounds(xx, yy, 0, 0));
    }
    gettimeofday(&end_time, nullptr); // record end time
    double buildTime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1e6; // 计算时间间隔
    std::cout << "Building Time :" << buildTime << "s" << std::endl;
//    expData << "Building Time :" << buildTime << "s" << std::endl;

    // calculate index size of R*tree
    x = tree.Query(RTree::AcceptAny(), Visitor());
    tree.indexSize(RTree::AcceptAny(), Visitor(), indexSize);
    resortCores(tree);

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

    cout << "Index Size of R-tree and Skyline Cores: " << indexSize / 1024.0 / 1024.0 << "MB" << endl;

    // calculate connectivity
    gettimeofday(&start_time, nullptr); // record start time
    connectivity(g, tree);
    gettimeofday(&end_time, nullptr); // record end time
    double connectTime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1e6; // 计算时间间隔
    std::cout << "Connecting Time :" << connectTime << "s" << std::endl;
//    expData << "Connecting Time :" << connectTime << "s" << std::endl;

    cout << "Index Size of R-tree and Skyline Cores: " << indexSize / 1024.0 / 1024.0 << "MB" << endl;
    int subgraphNum = 0;
    int connectEdgeNum = 0;
    for (auto &core: skylineCores) {
        subgraphNum += core.subGraphs.size();
        connectEdgeNum += core.connectEdges.size();
    }
    indexSize += height * subgraphNum * 4;
    indexSize += (outEdgeNum + connectEdgeNum) * 8;
    cout << "The number of skyline cores : " << skylineCores.size() << endl;

    cout << "Index Size of R-tree and Skyline Cores and Tags and Edges: " << indexSize / 1024.0 / 1024.0 << "MB"
         << endl;

    queryTest(tree, g);
//    expData << "Index Size of R-tree and Skyline Cores and Tags and Edges: " << indexSize / 1024.0 / 1024.0 << "MB"
//            << endl;

    return 0;
}