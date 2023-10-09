#include <string>
#include <ctime>
#include <sys/time.h>
#include <fstream>
#include <vector>
#include <map>
#include <stdio.h>
#include <queue>
#include <deque>
#include <omp.h>
#include <boost/multi_array.hpp>
#include "RStarTree.h"
#include "utilities/UnionFind.h"
#include "multilayer/multilayer_graph.h"

using namespace std;
const int dimensions = 9;
const int min_child_items = 32;
const int max_child_items = 64;
typedef RStarTree<int, dimensions, min_child_items, max_child_items> RTree;
typedef RTree::BoundingBox BoundingBox;

struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        std::hash<int> hasher;
        size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

int vertexNum;
int maxLayer;
vector<MultiSkylineCore> skylineCores;
map<vector<int>, int> mapCoreId;
vector<int> maxVec;
unordered_map<vector<int>, int, VectorHash> coreMax;
int coreID;
//boost::multi_array<int, dimensions> coreMatrix;
//boost::multi_array<int, dimensions> nodeMatrix;
unordered_map<vector<int>, int, VectorHash> mapCoCore;
unordered_map<vector<int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>, VectorHash> mapCoNode;
//map<vector<int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>> mapCoNode;
vector<vector<int>> traverseVec;
string dataset;
vector<vector<int>> possibleCores;

ofstream expData;

//BoundingBox bounds(int x, int y, int w, int h) {
//    BoundingBox bb;
//
//    bb.edges[0].first = x;
//    bb.edges[0].second = x + w;
//
//    bb.edges[1].first = y;
//    bb.edges[1].second = y + h;
//
//    return bb;
//}

BoundingBox MultiBounds(vector<int> vec) {
    BoundingBox bb;
    for (int i = 0; i < vec.size(); i++) {
        bb.edges[i].first = vec[i];
        bb.edges[i].second = vec[i];
    }
    return bb;
}

BoundingBox MultiBoundsQuery(vector<int> &vec) {
    BoundingBox bb;
    for (int i = 0; i < vec.size(); i++) {
        bb.edges[i].first = vec[i];
        bb.edges[i].second = vec[i] + maxVec[i];
    }
    return bb;
}

BoundingBox MultiBoundsConnect(vector<int> &vec) {
    BoundingBox bb;
    for (int i = 0; i < vec.size(); i++) {
        if (i != coreID) {
            bb.edges[i].first = vec[i];
            bb.edges[i].second = vec[i] + maxVec[i];
        }
    }
    bb.edges[coreID].first = vec[coreID];
    bb.edges[coreID].second = vec[coreID];
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
                   MultilayerGraph &multilayer_graph,  vector<bool> &isResult, UnionFind &UF) {
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
//        for (int i = 0; i < )
//    }

//    map<int, vector<int>> subGraph;
//    for (auto &vertex : results) {
//        subGraph[UF.Find(vertex)].push_back(vertex);
//    }
//    cout << " Connect Num:" << subGraph.size() << " ";
//    cout << " " << results.size() << " ";
}


void queryTest(RTree &tree, MultilayerGraph &multilayer_graph) {
    int test_case = 10000;
    vector<int> max_vec;
    max_vec.resize(maxLayer);
    clock_t start_time, end_time;

    int cnt = 0;
    double total_time = 0;
    for (int i = 0; i < test_case; i++) {
        vector<int> query_vec = possibleCores[rand() % possibleCores.size()];
//        query_vec = {2, 0, 0, 1, 4, 0, 0};
        cout << i << "/" << test_case;

        for (int j = 0; j < query_vec.size(); j++) {
            cout << " " << query_vec[j];
        }
//        query_vec.resize(maxLayer);
//        int vecMax = 10;
//        if (dataset == "homo") {
//            vecMax = 3;
//        } else if (dataset == "sacchcere") {
//            vecMax = 3;
//        } else if (dataset == "higgs") {
//            vecMax = 10;
//        } else if (dataset == "dblp") {
//            vecMax = 3;
//        }
//        cout << "Query vector:";
//        for (int j = 0; j < maxLayer; j++) {
//            query_vec[j] = rand() % vecMax;
////            query_vec[j] = 1;
//            cout << " " << query_vec[j];
//        }
//        if (mapCoreId.find(query_vec) == mapCoreId.end()) {
//            cout << endl;
//            continue;
//        }
        BoundingBox bound = MultiBoundsQuery(query_vec);
        vector<int> cores;
        vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *> nodes;
        double queryTime = 0.0;
        start_time = clock();
        tree.QueryPrune(RTree::AcceptEnclosing(bound), Visitor(), cores, nodes);
        end_time = clock();
        queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
//        total_time += queryTime;

        if (cores.size() > 0 || nodes.size() > 0) {
            vector<bool> isResult;
            isResult.resize(vertexNum);
            std::fill_n(isResult.begin(), isResult.size(), false);
            UnionFind UF(vertexNum);
            start_time = clock();
            getSubGraphUF(cores, nodes, multilayer_graph, isResult, UF);
            end_time = clock();
            queryTime += (double) (end_time - start_time) / CLOCKS_PER_SEC * 1000;
//            cout << queryTime << endl;
            cnt++;
        } else {
            cout << endl;
        }

        cnt++;
        cout << " " << queryTime << "ms" << endl;
        total_time += queryTime;
//        total_time += (double)queryTime / CLOCKS_PER_SEC * 1000;
    }
    cout << total_time / (double)cnt << "ms" << endl;
    expData << total_time / (double)cnt << "ms" << endl;
}


bool isEqual(std::vector<int> const& a, std::vector<int> const& b) {
    if (a.size() != b.size()) return false;
    return std::equal(a.begin(), a.end(), b.begin());
}

void preProcess(vector<MultiSkylineCore> &skylineCores, MultilayerGraph &multilayer_graph) {
    maxVec.resize(multilayer_graph.number_of_layers);
    int cnt = 0;
    for (int i = 0; i < skylineCores.size(); i++) {
        mapCoreId.insert({skylineCores[i].parameters, i});
        for (int j = 0; j < maxLayer; j++) {
            maxVec[j] = max(maxVec[j], skylineCores[i].parameters[j]);
        }
        printf("%d / %d\n", i, skylineCores.size());
    }

    multilayer_graph.neighbors.resize(vertexNum);
    for (int i = 0; i < multilayer_graph.number_of_nodes; i++) {
        printf("%d / %d\n", i, multilayer_graph.number_of_nodes);
        int vertexID = i + 1;
        set<int> neighbors;
        for (auto layer : multilayer_graph.adjacency_list[vertexID]) {
            for (auto neighbor : layer) {
                neighbors.insert(neighbor);
            }
        }
        int setSize = neighbors.size();
        multilayer_graph.neighbors[vertexID].reserve(setSize);
        for (auto it = neighbors.begin(); it != neighbors.end(); it++) {
            multilayer_graph.neighbors[vertexID].push_back(*it);
        }
    }
}

bool compareVec(const vector<int> &a, const vector<int> &b) {
    for (int i = 0; i < a.size(); i++) {
        if (a[i] != b[i]) {
            return a[i] > b[i];
        }
    }
}

void loadCoreMax(ifstream &ifstream1) {
    ifstream1 >> coreID;
    string line;
    while (getline(ifstream1, line)) {
        istringstream iss(line);
        vector<int> numbers;
        int number;
        int cnt = 0;
        while (iss >> number) {
            cnt++;
            if (cnt < dimensions)
                numbers.push_back(number);
        }
        if (numbers.size() == 0)
            continue;
        coreMax[numbers] = number;
        traverseVec.push_back(numbers);
        auto it = numbers.begin();
        it += coreID;
        numbers.insert(it, number);
        for (int i = 0; i <= number; i++) {
            numbers[coreID] = i;
            possibleCores.push_back(numbers);
        }
//        numbers.push_back(number);
//        swap(numbers[coreID], numbers[numbers.size()-1]);
    }
    std::sort(traverseVec.begin(), traverseVec.end(), compareVec);
//    // initialize matrix
//    for (int i = 0; i < dimensions; i++) {
//        if (i == 0) {
//            int tempMax = 0;
//            for (auto &item : coresPara) {
//                tempMax = max(item[0], tempMax);
//            }
//            coreMatrix.resize(tempMax + 1);
//            nodeMatrix.resize(tempMax + 1);
//        } else {
//            vector<int> tempInt;
//            tempInt.resize(100);
//            nodeMatrix[0] = tempInt;
//        }
//    }
    cout << "x";
}

void loadSkylineCore(string skylinePath) {
    ifstream inFile(skylinePath);
    int skylineNum;
    inFile >> skylineNum;
    skylineCores.resize(skylineNum);
    for (int i = 0; i < skylineNum; i++) {
        MultiSkylineCore tempCore;
        vector<int> tempVec;
        tempVec.resize(maxLayer);
        for (int j = 0; j < maxLayer; j++) {
            int tempNum;
            inFile >> tempNum;
            tempVec[j] = tempNum;
        }
        tempCore.parameters = tempVec;
        int tempVertexNum;
        inFile >> tempVertexNum;
        for (int j = 0; j < tempVertexNum; j++) {
            int tempID;
            inFile >> tempID;
            tempCore.vertices.push_back(tempID);
        }
        skylineCores[i] = tempCore;
    }
    inFile.close();
}

void rtreeDFS(RStarTree<int, dimensions, min_child_items, max_child_items>::Node * node, int &coreCnt, vector<MultiSkylineCore> &newSkylineCores) {
    vector<int> coordinate;
    for (int i = 0; i < dimensions; i++) {
        coordinate.push_back(node->bound.edges[i].first);
    }
    mapCoNode[coordinate].push_back(node);
//    nodeMatrix[node->bound.edges[0].first][node->bound.edges[1].first] = 1;
    if (node->hasLeaves) {
        node->left = coreCnt;
        for (auto item1 : node->items) {
            vector<int> tempVec;
            for (int i = 0; i < dimensions; i++) {
                tempVec.push_back(item1->bound.edges[i].first);
            }
            static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Leaf *>(item1)->coreId = coreCnt;
            newSkylineCores.push_back(skylineCores[mapCoreId[tempVec]]);
            newSkylineCores.back().id = coreCnt;
//            mapCoCore[tempVec] = coreCnt;
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
    vector<MultiSkylineCore> newSkylineCores;
//    map<pair<int, int>, int> newMapCoreId;
    if (node != NULL) {
        rtreeDFS(node, coreCnt, newSkylineCores);
    }
    skylineCores = newSkylineCores;
    for (int i = 0; i < skylineCores.size(); i++) {
        mapCoCore[skylineCores[i].parameters] = i;
    }
}

bool compareSkylineCore(const MultiSkylineCore &a, const MultiSkylineCore &b) {
    if (a.parameters.back() != b.parameters.back()) {
        return a.parameters.back() > b.parameters.back();
    } else {
        return a.parameters[0] < b.parameters[0];
    }
}

int outEdgeNum = 0;

void connectivity2(MultilayerGraph &multilayer_graph, RTree &tree) {
    timeval start_time, end_time;
    int cnt = 0;
    UnionFind UF(vertexNum);
    UnionFind UF2(vertexNum);
    vector<bool> nodes(vertexNum);
    std::fill_n(nodes.begin(), nodes.size(), false);
    vector<int> preResult;
    int traverseCnt = 0;
    for (int i = 0; i < traverseVec.size(); i++) {
        cout << i << "/" << traverseVec.size() << endl;
        int nowMax = coreMax[traverseVec[i]];
        vector<int> nowPoint = traverseVec[i];
        nowPoint.insert(nowPoint.begin() + coreID, nowMax);
        if (i != 0 && coreMax[traverseVec[i]] == coreMax[traverseVec[i-1]]) {
            bool isBlank = true;
            for (int j = 0; j <= nowMax; j++) {
                nowPoint[coreID] = j;
                if (mapCoCore.find(nowPoint) != mapCoCore.end() || mapCoNode.find(nowPoint) != mapCoNode.end()) {
                    isBlank = false;
                    break;
                }
            }
            if (isBlank) {
                continue;
            }
        }
        UF.Init(preResult);
        UF2.Init(preResult);
        for (auto &element : preResult) {
            nodes[element] = false;
        }
        preResult.clear();
        for (int tempNum = nowMax; tempNum >= 0; tempNum--) {
            nowPoint[coreID] = tempNum;
//            bool flag = true;
//            for (auto element : nowPoint) {
//                if (element != 0) {
//                    flag = false;
//                    break;
//                }
//            }
//            if (flag) {
//                break;
//            }
            vector<int> cores;
            vector<int> tempResults;
            // get new cores
            BoundingBox bound;
            bound = MultiBoundsConnect(nowPoint);
            tree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);

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
                for (int node: multilayer_graph.neighbors[tempResults[j]]) {
                    if (nodes[node] && UF.Find(tempResults[j]) != UF.Find(node)) {
                        UF.Unite(tempResults[j], node);
                    }
                }
            }

            //// Count connected subgraphs of current core
            if (mapCoCore.find(nowPoint) != mapCoCore.end()) {
                int tempCoreID = mapCoCore[nowPoint];
                map<int, vector<int>> rootVec;
                for (auto vertex: skylineCores[tempCoreID].vertices) {
                    int tempRoot = UF.Find(vertex);
                    rootVec[tempRoot].push_back(vertex);
                }
                for (auto &item1: rootVec) {
                    skylineCores[tempCoreID].subGraphs.push_back(item1.second);
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
                for (int node: multilayer_graph.neighbors[tempResults[j]]) {
                    if (nodes[node] && UF2.Find(tempResults[j]) != UF2.Find(node)) {
                        UF2.Unite(tempResults[j], node);
                        int tempCoreId = newVertexCore[tempResults[j]];
                        skylineCores[tempCoreId].connectEdges.push_back({tempResults[j], node});
                    }
                }
            }

            //// calculate tags and out edges for R-tree node
            if (mapCoNode.find(nowPoint) != mapCoNode.end()) {
//                cout << "!!";
                auto RNodes = mapCoNode.find(nowPoint)->second;
                vector<bool> tempNodes;
                tempNodes.resize(vertexNum);
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

//  dataset = "higgs";
//    dataset = "homo";
//  dataset = "sacchcere";
//  dataset = "dblp";
//  dataset = "example";
    dataset = argv[1];
    string expPath = "./multilayer/" + dataset + "/" + dataset + "-exp.txt";
    expData.open(expPath);
    string datasetPath = "./multilayer/" + dataset + "/" + dataset + ".txt";
    //// load graph
    cout << "Loading graph...." << endl;
    MultilayerGraph multilayer_graph(datasetPath);
    maxLayer = multilayer_graph.number_of_layers;
    vertexNum = multilayer_graph.number_of_nodes;
    vertexNum++;

    //// load skyline cores
    ifstream inFile;
    string skylinePath = "./multilayer/" + dataset + "/" + dataset + "_skyline";
    cout << "Loading skyline cores..." << endl;
    loadSkylineCore(skylinePath);

    //// preprocess
    preProcess(skylineCores, multilayer_graph);

    cout << skylineCores.size() << endl;

    string filePath = "./multilayer/" + dataset + "/" + dataset + "_numMax";
    inFile.open(filePath);
    loadCoreMax(inFile);
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
        tree.Insert(i, MultiBounds(skylineCores[i].parameters));
    }
    gettimeofday(&end_time, nullptr); // record end time
    double buildTime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1e6; // 计算时间间隔
    std::cout << "Building Time :" << buildTime << "s" << std::endl;
    expData << "Building Time :" << buildTime << "s" << std::endl;

    // calculate index size of R*tree
    tree.indexSizeMulti(RTree::AcceptAny(), Visitor(), indexSize);
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

    gettimeofday(&start_time, nullptr); // record start time
//    connectivity(multilayer_graph, tree);
    connectivity2(multilayer_graph, tree);
    gettimeofday(&end_time, nullptr); // record end time
    double connectTime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1e6; // 计算时间间隔
    std::cout << "Connect Time :" << connectTime << "s" << std::endl;
    expData << "Connect Time :" << connectTime << "s" << std::endl;

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

//    query
    queryTest(tree, multilayer_graph);
//    queryTestWithPoint(tree, g);
    expData.close();
    return 0;
}
