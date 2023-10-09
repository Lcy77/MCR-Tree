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
const int dimensions = 4;
const int min_child_items = 32;
const int max_child_items = 64;
typedef RStarTree<int, dimensions, min_child_items, max_child_items> RTree;
typedef RTree::BoundingBox BoundingBox;

struct VectorHash {
    size_t operator()(const std::vector<int> &v) const {
        std::hash<int> hasher;
        size_t seed = 0;
        for (int i: v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
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

struct UpdatedVertex {
    int vertexID;
    vector<vector<int>> newCores;
    vector<vector<int>> oldCores;
};

struct addEdge {
    vector<UpdatedVertex> updatedVertices;
    vector<vector<int>> vertexPairs;  // from vertex; to vertex; layer
};

struct deleteEdge {
    vector<UpdatedVertex> updatedVertices;
    vector<vector<int>> vertexPairs;
};


bool isEqual(std::vector<int> const &a, std::vector<int> const &b) {
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
        for (auto layer: multilayer_graph.adjacency_list[vertexID]) {
            for (auto neighbor: layer) {
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

void rtreeDFS(RStarTree<int, dimensions, min_child_items, max_child_items>::Node *node, int &coreCnt,
              vector<MultiSkylineCore> &newSkylineCores) {
    vector<int> coordinate;
    for (int i = 0; i < dimensions; i++) {
        coordinate.push_back(node->bound.edges[i].first);
    }
    mapCoNode[coordinate].push_back(node);
//    nodeMatrix[node->bound.edges[0].first][node->bound.edges[1].first] = 1;
    if (node->hasLeaves) {
        node->left = coreCnt;
        for (auto item1: node->items) {
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

void updatedRtreeDFS(RStarTree<int, dimensions, min_child_items, max_child_items>::Node *node, int &coreCnt,
                     unordered_map<vector<int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>, VectorHash> &updatedMapCoNode,
                     unordered_map<vector<int>, int, VectorHash> &updatedMapCoCore,
                     unordered_map<vector<int>, int, VectorHash> &newMapCoCore,
                     vector<MultiSkylineCore> &updatedSkylineCores,
                     vector<MultiSkylineCore> &newSkylineCores) {
    vector<int> coordinate;
    for (int i = 0; i < dimensions; i++) {
        coordinate.push_back(node->bound.edges[i].first);
    }
    updatedMapCoNode[coordinate].push_back(node);
    if (node->hasLeaves) {
        node->left = coreCnt;
        for (auto item1: node->items) {
            vector<int> tempVec;
            for (int i = 0; i < dimensions; i++) {
                tempVec.push_back(item1->bound.edges[i].first);
            }
            static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Leaf *>(item1)->coreId = coreCnt;
            newSkylineCores.push_back(updatedSkylineCores[updatedMapCoCore[tempVec]]);
            newSkylineCores.back().id = coreCnt;
            newMapCoCore[tempVec] = coreCnt;
            coreCnt++;
        }
        node->right = coreCnt - 1;
    } else {
        for (auto item1: node->items) {
            updatedRtreeDFS(static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(item1),
                            coreCnt, updatedMapCoNode, updatedMapCoCore, newMapCoCore,
                            updatedSkylineCores, newSkylineCores);
        }
        node->left = static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(node->items[0])->left;
        node->right = static_cast<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>(node->items.back())->right;
//        cout << node->left << "!!";
    }
//    cout << "range1:" << node->left << " " << node->right << endl;

}

void updatedResortCores(RTree &tree,
                        unordered_map<vector<int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>, VectorHash> &updatedMapCoNode,
                        unordered_map<vector<int>, int, VectorHash> &updatedMapCoCore,
                        vector<MultiSkylineCore> &updatedSkylineCores,
                        map<int, vector<vector<int>>> &updatedMapCoreNode) {
    auto node = tree.m_root;
    int coreCnt = 0;
    unordered_map<vector<int>, int, VectorHash> newMapCoCore;
    vector<MultiSkylineCore> newSkylineCores;
    if (node != NULL) {
        updatedRtreeDFS(node, coreCnt, updatedMapCoNode, updatedMapCoCore, newMapCoCore,
                        updatedSkylineCores, newSkylineCores);
    }
    updatedSkylineCores = newSkylineCores;
    updatedMapCoCore = newMapCoCore;
    for (auto item: updatedMapCoNode) {
        for (auto treeNode: item.second) {
            vector<int> tempVec;
            for (int i = 0; i < dimensions; i++) {
                tempVec.push_back(treeNode->bound.edges[i].first);
            }
            for (int i = treeNode->left; i <= treeNode->right; i++) {
                updatedMapCoreNode[i].push_back(tempVec);
            }
        }
    }
}

bool compareSkylineCore(const MultiSkylineCore &a, const MultiSkylineCore &b) {
    if (a.parameters.back() != b.parameters.back()) {
        return a.parameters.back() > b.parameters.back();
    } else {
        return a.parameters[0] < b.parameters[0];
    }
}

void addOne(RTree &tree, MultilayerGraph &multilayer_graph, int num) {
    string filePath = "./multilayer/" + dataset + "/" + dataset + "_add" + to_string(num);
    ifstream updateFile(filePath);
    string line;
    vector<addEdge> testCase;
    while (std::getline(updateFile, line)) {
        if (line.size() <= 0)
            break;
        if (line[0] == 'A') {
            int pos = line.find(":");
            line = line.substr(pos + 2);
            istringstream iss(line);
            int fromVertex, toVertex;
            int layer;
            addEdge tempAdd;
            while (iss >> fromVertex >> toVertex >> layer) {
//                cout << fromVertex << " " << toVertex << endl;
                vector<int> tempPairs(3);
                tempPairs[0] = fromVertex;
                tempPairs[1] = toVertex;
                tempPairs[2] = layer;
                tempAdd.vertexPairs.push_back(tempPairs);
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
                // Extract the core's num
                stringstream ss(token);
                vector<int> tempCore;
                string tempNum;
                while (getline(ss, tempNum, ',')) {
                    tempCore.push_back(stoi(tempNum));
                }
                tempVertex.newCores.push_back(tempCore);
                strCoordinatesNew.erase(0, pos + delimiter.length()); // Remove the processed part
            }
            stringstream ss1(strCoordinatesNew);
            vector<int> tempCore1;
            string tempNum1;
            while (getline(ss1, tempNum1, ',')) {
                tempCore1.push_back(stoi(tempNum1));
            }
            tempVertex.newCores.push_back(tempCore1);

            //// Process the old coordinates
            strCoordinatesOld = strCoordinatesOld.substr(1, strCoordinatesOld.length() - 2);
            pos = 0;
            while ((pos = strCoordinatesOld.find(delimiter)) != std::string::npos) {
                token = strCoordinatesOld.substr(0, pos); // Extract a coordinate pair string
                // Extract the core's num
                stringstream ss(token);
                vector<int> tempCore;
                string tempNum;
                while (getline(ss, tempNum, ',')) {
                    tempCore.push_back(stoi(tempNum));
                }
                tempVertex.oldCores.push_back(tempCore);
                strCoordinatesOld.erase(0, pos + delimiter.length()); // Remove the processed part
            }
            stringstream ss2(strCoordinatesOld);
            vector<int> tempCore2;
            string tempNum2;
            while (getline(ss2, tempNum2, ',')) {
                tempCore2.push_back(stoi(tempNum2));
            }
            tempVertex.oldCores.push_back(tempCore2);
            testCase[testCase.size() - 1].updatedVertices.push_back(tempVertex);
        }
    }
    updateFile.close();

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
        auto updatedSkylineCores = skylineCores;
        unordered_map<vector<int>, int, VectorHash> updatedMapCoCore = mapCoCore;
        unordered_map<vector<int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>, VectorHash> updatedMapCoNode;
        RTree originTree;
        for (int i = 0; i < updatedSkylineCores.size(); i++) {
            originTree.Insert(i, MultiBounds(updatedSkylineCores[i].parameters));
        }

        map<vector<int>, vector<int>> mapOldSkylineCore;
        map<vector<int>, vector<int>> mapNewSkylineCore;
        for (auto &tempVertex: updatedVertices) {
            for (int i = 0; i < tempVertex.oldCores.size(); i++) {
                mapOldSkylineCore[tempVertex.oldCores[i]].push_back(tempVertex.vertexID);
            }
            for (int i = 0; i < tempVertex.newCores.size(); i++) {
                mapNewSkylineCore[tempVertex.newCores[i]].push_back(tempVertex.vertexID);
            }
        }
        // build 2 R-tree for updated vertices
        vector<MultiSkylineCore> oldSkylineCore;
        vector<MultiSkylineCore> newSkylineCore;
        RTree oldTree;
        RTree newTree;
        map<vector<int>, int> mapOldCoreId;
        map<vector<int>, int> mapNewCoreId;
        int index1 = 0;
        int index2 = 0;
        for (auto &core: mapOldSkylineCore) {
            MultiSkylineCore tempSkylineCore = MultiSkylineCore();
            oldTree.Insert(index1, MultiBounds(core.first));
            tempSkylineCore.parameters = core.first;
            tempSkylineCore.vertices = core.second;
            std::sort(tempSkylineCore.vertices.begin(), tempSkylineCore.vertices.end());
            oldSkylineCore.push_back(tempSkylineCore);
            mapOldCoreId[tempSkylineCore.parameters] = index1;
            index1++;
        }
        for (auto &core: mapNewSkylineCore) {
            MultiSkylineCore tempSkylineCore = MultiSkylineCore();
            newTree.Insert(index2, MultiBounds(core.first));
            tempSkylineCore.parameters = core.first;
            tempSkylineCore.vertices = core.second;
            std::sort(tempSkylineCore.vertices.begin(), tempSkylineCore.vertices.end());
            newSkylineCore.push_back(tempSkylineCore);
            mapNewCoreId[tempSkylineCore.parameters] = index2;
            index2++;
        }
        long long size = 0;
        oldTree.indexSizeMulti(RTree::AcceptAny(), Visitor(), size, mapOldCoreId);
        newTree.indexSizeMulti(RTree::AcceptAny(), Visitor(), size, mapNewCoreId);

        // updated original skyline core
        for (auto vertex: updatedVertices) {
            // find the changed skyline cores of vertex
            std::sort(vertex.newCores.begin(), vertex.newCores.end());
            std::sort(vertex.oldCores.begin(), vertex.oldCores.end());
            vector<vector<int>> coresAdded;
            vector<vector<int>> coresDeleted;
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
                    BoundingBox bound = MultiBounds(coresDeleted[i]);
                    originTree.RemoveBoundedArea(bound);
                    updatedMapCoCore.erase(coresDeleted[i]);
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
                auto &core = coresAdded[i];
                if (updatedMapCoCore.find(core) == updatedMapCoCore.end()) {
                    MultiSkylineCore tempSkylineCore = MultiSkylineCore();
                    tempSkylineCore.parameters = core;
                    updatedSkylineCores.push_back(tempSkylineCore);
                    updatedMapCoCore[core] = updatedSkylineCores.size() - 1;
                    BoundingBox bound = MultiBounds(core);
                    originTree.Insert(updatedSkylineCores.size() - 1, bound);
                }
                int coreId = updatedMapCoCore[core];
                updatedSkylineCores[coreId].vertices.push_back(vertex.vertexID);
                updatedSkylineCores[coreId].addedVertices.push_back(vertex.vertexID);
            }
        }
        map<int, vector<vector<int>>> mapCoreNode;  // coreID -> nodes
        updatedResortCores(originTree, updatedMapCoNode, updatedMapCoCore, updatedSkylineCores, mapCoreNode);

        int updatedTraverseStart = 0;
        vector<int> updatedTraverseVertexNum(traverseVec.size() + 1);
        set<int> tempUpdatedVertex;
        bool flagStart = false;
        for (int i = 0; i < traverseVec.size(); i++) {
            vector<int> nowPoint = traverseVec[i];
            int nowMax = coreMax[traverseVec[i]];
            nowPoint.insert(nowPoint.begin() + coreID, nowMax);

            for (int tempNum = nowMax; tempNum >= 0; tempNum--) {
                nowPoint[coreID] = tempNum;
                vector<int> tempOldCores;
                vector<int> tempNewCores;
                oldTree.QueryCore(RTree::AcceptEnclosing(MultiBounds(nowPoint)), Visitor(), tempOldCores);
                newTree.QueryCore(RTree::AcceptEnclosing(MultiBounds(nowPoint)), Visitor(), tempNewCores);
                if (!flagStart && tempOldCores.size() != 0) {
                    for (auto &vertex: oldSkylineCore[tempOldCores[0]].vertices) {
                        tempUpdatedVertex.insert(vertex);
                    }
                }
                if (!flagStart && tempNewCores.size() != 0) {
                    if (tempOldCores.size() != 0) {
                        if (oldSkylineCore[tempOldCores[0]].vertices != newSkylineCore[tempNewCores[0]].vertices) {
                            updatedTraverseStart = i;
                            flagStart = true;
                        }
                    } else {
                        updatedTraverseStart = i;
                        flagStart = true;
                    }
                }
            }
            updatedTraverseVertexNum[i] = tempUpdatedVertex.size();
        }

//        cout << "updatedTraverseStart" << updatedTraverseStart << endl;
        start_time = clock();
        set<int> changedCoreId;
        set<vector<int>> changedNodeCo;

        UnionFind UF(vertexNum);
        vector<bool> nodes(vertexNum);
        std::fill_n(nodes.begin(), nodes.size(), false);
        vector<int> preResult;
        vector<bool> oldFlag(vertexNum);
        vector<int> preOld;
        for (int i = updatedTraverseStart; i < traverseVec.size(); i++) {
            int nowMax = coreMax[traverseVec[i]];
            vector<int> nowPoint = traverseVec[i];
            nowPoint.insert(nowPoint.begin() + coreID, nowMax);
            if (i != 0 && coreMax[traverseVec[i]] == coreMax[traverseVec[i - 1]]) {
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
            for (auto &element: preResult) {
                nodes[element] = false;
            }
            preResult.clear();
            BoundingBox bound;
            bool start = false;
            int oldCnt = 0;
            set<int> addedVertices;
            for (auto &element: preOld) {
                oldFlag[element] = false;
            }
            preOld.clear();
            for (int tempNum = nowMax; tempNum >= 0; tempNum--) {
                vector<int> oldCores;
                vector<int> newCores;
                vector<int> tempResults;
                nowPoint[coreID] = tempNum;
                BoundingBox tempBound = MultiBoundsConnect(nowPoint);;
                oldTree.QueryCore(RTree::AcceptEnclosing(tempBound), Visitor(), oldCores);
                newTree.QueryCore(RTree::AcceptEnclosing(tempBound), Visitor(), newCores);
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
                if (addedVertices.size() != 0) {
                    //// reconstruct the subgraphs
                    vector<int> cores;
                    if (!start) {
                        bound = MultiBoundsQuery(nowPoint);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        bound = MultiBoundsConnect(nowPoint);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
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

                    // Check whether the new point can change the connectivity
                    vector<pair<int, int>> extraEdges;
                    auto it = addedVertices.begin();
                    while (it != addedVertices.end()) {
                        int vertexId = *it;
                        int root1 = UF.Find(vertexId);
                        for (int node: multilayer_graph.neighbors[vertexId]) {
                            if (nodes[node]) {
                                int root2 = UF.Find(node);
                                if (root1 != root2 && addedVertices.find(node) == addedVertices.end()) {
                                    extraEdges.push_back({vertexId, node});
                                }
                                UF.Unite(vertexId, node);
                            }
                        }
                        it++;
                    }

                    if (updatedMapCoCore.find(nowPoint) != updatedMapCoCore.end()) {
                        // updated the current core
                        bool changed = false;
                        int tempCoreID = updatedMapCoCore[nowPoint];
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
                        bound = MultiBoundsQuery(nowPoint);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        bound = MultiBoundsConnect(nowPoint);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
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
                        if (UF.Find(vertexPair[0]) != UF.Find(vertexPair[1])) {
                            change = true;
                            UF.Unite(vertexPair[0], vertexPair[1]);
                            updatedSkylineCores[newVertexCore[vertexPair[0]]].connectEdges.push_back(
                                    {vertexPair[0], vertexPair[1]});
                            for (auto nodeCore: mapCoreNode[newVertexCore[vertexPair[0]]]) {
                                changedNodeCo.insert(nodeCore);
                            }
                        }
                    }
//                    if (oldCnt >= updatedKVertexNum[tempK]) {
//                        break;
//                    }
                    if (!change && oldCnt >= updatedTraverseVertexNum[i]) {
                        break;
                    }
                }
                // updated R-tree node
                if (changedNodeCo.find(nowPoint) != changedNodeCo.end()) {
//                    cout << "!!";
                    auto RNodes = updatedMapCoNode.find(nowPoint)->second;
                    vector<bool> tempNodes;
                    tempNodes.resize(vertexNum);
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
                    changedNodeCo.erase(nowPoint);
                }
            }
        }
        end_time = clock();
        total_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
    }
    cout << "Update time : " << total_time / (double) testCase.size() << "s" << endl;
}

void deleteOne(RTree &tree, MultilayerGraph &multilayer_graph, int num) {
    string filePath = "./multilayer/" + dataset + "/" + dataset + "_delete" + to_string(num);
    ifstream updateFile(filePath);
    string line;
    vector<addEdge> testCase;
    while (std::getline(updateFile, line)) {
        if (line.size() <= 0)
            break;
        if (line[0] == 'D') {
            int pos = line.find(":");
            line = line.substr(pos + 2);
            istringstream iss(line);
            int fromVertex, toVertex;
            int layer;
            addEdge tempAdd;
            while (iss >> fromVertex >> toVertex >> layer) {
//                cout << fromVertex << " " << toVertex << endl;
                vector<int> tempPairs(3);
                tempPairs[0] = fromVertex;
                tempPairs[1] = toVertex;
                tempPairs[2] = layer;
                tempAdd.vertexPairs.push_back(tempPairs);
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
                // Extract the core's num
                stringstream ss(token);
                vector<int> tempCore;
                string tempNum;
                while (getline(ss, tempNum, ',')) {
                    tempCore.push_back(stoi(tempNum));
                }
                tempVertex.newCores.push_back(tempCore);
                strCoordinatesNew.erase(0, pos + delimiter.length()); // Remove the processed part
            }
            stringstream ss1(strCoordinatesNew);
            vector<int> tempCore1;
            string tempNum1;
            while (getline(ss1, tempNum1, ',')) {
                tempCore1.push_back(stoi(tempNum1));
            }
            tempVertex.newCores.push_back(tempCore1);

            //// Process the old coordinates
            strCoordinatesOld = strCoordinatesOld.substr(1, strCoordinatesOld.length() - 2);
            pos = 0;
            while ((pos = strCoordinatesOld.find(delimiter)) != std::string::npos) {
                token = strCoordinatesOld.substr(0, pos); // Extract a coordinate pair string
                // Extract the core's num
                stringstream ss(token);
                vector<int> tempCore;
                string tempNum;
                while (getline(ss, tempNum, ',')) {
                    tempCore.push_back(stoi(tempNum));
                }
                tempVertex.oldCores.push_back(tempCore);
                strCoordinatesOld.erase(0, pos + delimiter.length()); // Remove the processed part
            }
            stringstream ss2(strCoordinatesOld);
            vector<int> tempCore2;
            string tempNum2;
            while (getline(ss2, tempNum2, ',')) {
                tempCore2.push_back(stoi(tempNum2));
            }
            tempVertex.oldCores.push_back(tempCore2);
            testCase[testCase.size() - 1].updatedVertices.push_back(tempVertex);
        }
    }
    updateFile.close();

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
        auto updatedSkylineCores = skylineCores;
        unordered_map<vector<int>, int, VectorHash> updatedMapCoCore = mapCoCore;
        unordered_map<vector<int>, vector<RStarTree<int, dimensions, min_child_items, max_child_items>::Node *>, VectorHash> updatedMapCoNode;
        RTree originTree;
        for (int i = 0; i < updatedSkylineCores.size(); i++) {
            originTree.Insert(i, MultiBounds(updatedSkylineCores[i].parameters));
        }

        map<vector<int>, vector<int>> mapOldSkylineCore;
        map<vector<int>, vector<int>> mapNewSkylineCore;
        for (auto &tempVertex: updatedVertices) {
            for (int i = 0; i < tempVertex.oldCores.size(); i++) {
                mapOldSkylineCore[tempVertex.oldCores[i]].push_back(tempVertex.vertexID);
            }
            for (int i = 0; i < tempVertex.newCores.size(); i++) {
                mapNewSkylineCore[tempVertex.newCores[i]].push_back(tempVertex.vertexID);
            }
        }
        auto updatedGraph = multilayer_graph;
        for (auto deleteEdge: tempCase.vertexPairs) {
            auto it = std::find(updatedGraph.adjacency_list[deleteEdge[0]][deleteEdge[2]].begin(),
                                updatedGraph.adjacency_list[deleteEdge[0]][deleteEdge[2]].end(), deleteEdge[1]);
            updatedGraph.adjacency_list[deleteEdge[0]][deleteEdge[2]].erase(it);
            it = std::find(updatedGraph.adjacency_list[deleteEdge[1]][deleteEdge[2]].begin(),
                           updatedGraph.adjacency_list[deleteEdge[1]][deleteEdge[2]].end(), deleteEdge[0]);
            updatedGraph.adjacency_list[deleteEdge[1]][deleteEdge[2]].erase(it);
        }
        updatedGraph.neighbors.clear();
        updatedGraph.neighbors.resize(vertexNum);
        for (int i = 0; i < updatedGraph.number_of_nodes; i++) {
//            printf("%d / %d\n", i, updatedGraph.number_of_nodes);
            int vertexID = i;
            if (dataset != "higgs")
                vertexID++;
            set<int> neighbors;
            for (auto layer: updatedGraph.adjacency_list[vertexID]) {
                for (auto neighbor: layer) {
                    neighbors.insert(neighbor);
                }
            }
            int setSize = neighbors.size();
            updatedGraph.neighbors[vertexID].reserve(setSize);
            for (auto it = neighbors.begin(); it != neighbors.end(); it++) {
                updatedGraph.neighbors[vertexID].push_back(*it);
            }
        }
        // build 2 R-tree for updated vertices
        vector<MultiSkylineCore> oldSkylineCore;
        vector<MultiSkylineCore> newSkylineCore;
        RTree oldTree;
        RTree newTree;
        map<vector<int>, int> mapOldCoreId;
        map<vector<int>, int> mapNewCoreId;
        int index1 = 0;
        int index2 = 0;
        for (auto &core: mapOldSkylineCore) {
            MultiSkylineCore tempSkylineCore = MultiSkylineCore();
            oldTree.Insert(index1, MultiBounds(core.first));
            tempSkylineCore.parameters = core.first;
            tempSkylineCore.vertices = core.second;
            oldSkylineCore.push_back(tempSkylineCore);
            mapOldCoreId[tempSkylineCore.parameters] = index1;
            index1++;
        }
        for (auto &core: mapNewSkylineCore) {
            MultiSkylineCore tempSkylineCore = MultiSkylineCore();
            newTree.Insert(index2, MultiBounds(core.first));
            tempSkylineCore.parameters = core.first;
            tempSkylineCore.vertices = core.second;
            newSkylineCore.push_back(tempSkylineCore);
            mapNewCoreId[tempSkylineCore.parameters] = index2;
            index2++;
        }
        long long size = 0;
        oldTree.indexSizeMulti(RTree::AcceptAny(), Visitor(), size, mapOldCoreId);

        newTree.indexSizeMulti(RTree::AcceptAny(), Visitor(), size, mapNewCoreId);

        // updated original skyline core
        for (auto vertex: updatedVertices) {
            // find the changed skyline cores of vertex
            std::sort(vertex.newCores.begin(), vertex.newCores.end());
            std::sort(vertex.oldCores.begin(), vertex.oldCores.end());
            vector<vector<int>> coresAdded;
            vector<vector<int>> coresDeleted;
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
                    BoundingBox bound = MultiBounds(coresDeleted[i]);
                    originTree.RemoveBoundedArea(bound);
                    updatedMapCoCore.erase(coresDeleted[i]);
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
                auto &core = coresAdded[i];
                if (updatedMapCoCore.find(core) == updatedMapCoCore.end()) {
                    MultiSkylineCore tempSkylineCore = MultiSkylineCore();
                    tempSkylineCore.parameters = core;
                    updatedSkylineCores.push_back(tempSkylineCore);
                    updatedMapCoCore[core] = updatedSkylineCores.size() - 1;
                    BoundingBox bound = MultiBounds(core);
                    originTree.Insert(updatedSkylineCores.size() - 1, bound);
                }
                int coreId = updatedMapCoCore[core];
                updatedSkylineCores[coreId].vertices.push_back(vertex.vertexID);
                updatedSkylineCores[coreId].addedVertices.push_back(vertex.vertexID);
            }
        }
        map<int, vector<vector<int>>> mapCoreNode;  // coreID -> nodes
        updatedResortCores(originTree, updatedMapCoNode, updatedMapCoCore, updatedSkylineCores, mapCoreNode);

        int updatedTraverseStart = 0;
        vector<int> updatedTraverseVertexNum(traverseVec.size() + 1);
        set<int> tempUpdatedVertex;
        bool flagStart = false;
        for (int i = 0; i < traverseVec.size(); i++) {
            vector<int> nowPoint = traverseVec[i];
            int nowMax = coreMax[traverseVec[i]];
            nowPoint.insert(nowPoint.begin() + coreID, nowMax);
            for (int tempNum = nowMax; tempNum >= 0; tempNum--) {
                nowPoint[coreID] = tempNum;
                vector<int> tempOldCores;
                vector<int> tempNewCores;
                oldTree.QueryCore(RTree::AcceptEnclosing(MultiBounds(nowPoint)), Visitor(), tempOldCores);
                newTree.QueryCore(RTree::AcceptEnclosing(MultiBounds(nowPoint)), Visitor(), tempNewCores);
                if (tempOldCores.size() != 0) {
                    for (auto &vertex: oldSkylineCore[tempOldCores[0]].vertices) {
                        tempUpdatedVertex.insert(vertex);
                    }
                }
                if (!flagStart && tempNewCores.size() != 0) {
                    if (tempOldCores.size() != 0) {
                        if (oldSkylineCore[tempOldCores[0]].vertices != newSkylineCore[tempNewCores[0]].vertices) {
                            updatedTraverseStart = i;
                            flagStart = true;
                        }
                    } else {
                        updatedTraverseStart = i;
                        flagStart = true;
                    }
                }
            }
            updatedTraverseVertexNum[i] = tempUpdatedVertex.size();
        }

        start_time = clock();
        set<int> changedCoreId;
        set<vector<int>> changedNodeCo;

        UnionFind UF(vertexNum);
        vector<bool> nodes(vertexNum);
        std::fill_n(nodes.begin(), nodes.size(), false);
        vector<int> preResult;
        vector<bool> newFlag(vertexNum);
        vector<int> preNew;
        for (int i = updatedTraverseStart; i < traverseVec.size(); i++) {
            int nowMax = coreMax[traverseVec[i]];
            vector<int> nowPoint = traverseVec[i];
            nowPoint.insert(nowPoint.begin() + coreID, nowMax);
            if (i != 0 && coreMax[traverseVec[i]] == coreMax[traverseVec[i - 1]]) {
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
            for (auto &element: preResult) {
                nodes[element] = false;
            }
            preResult.clear();
            BoundingBox bound;
            bool start = false;
            int newCnt = 0;
            for (auto &element: preNew) {
                newFlag[element] = false;
            }
            preNew.clear();
            set<int> deletedVertices;
            for (int tempNum = nowMax; tempNum >= 0; tempNum--) {
                vector<int> oldCores;
                vector<int> newCores;
                vector<int> tempResults;
                nowPoint[coreID] = tempNum;
                BoundingBox tempBound = MultiBoundsConnect(nowPoint);;
                oldTree.QueryCore(RTree::AcceptEnclosing(tempBound), Visitor(), oldCores);
                newTree.QueryCore(RTree::AcceptEnclosing(tempBound), Visitor(), newCores);
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
                if (deletedVertices.size() != 0) {
                    //// reconstruct the subgraphs
                    vector<int> cores;
                    if (!start) {
                        bound = MultiBoundsQuery(nowPoint);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        bound = MultiBoundsConnect(nowPoint);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
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

                    //// Calculate connectivity
                    for (int j = 0; j < tempResults.size(); j++) {
                        for (int node: updatedGraph.neighbors[tempResults[j]]) {
                            if (nodes[node]) {
                                UF.Unite(tempResults[j], node);
                            }
                        }
                    }

                    bool change = false;
                    if (mapCoCore.find(nowPoint) != updatedMapCoCore.end()) {
                        // updated the current core
                        int coreId = updatedMapCoCore[nowPoint];
                        // Determine whether the deleted point can change the connectivity
                        for (int j = 0; j < updatedSkylineCores[coreId].subGraphs.size(); j++) {
                            map<int, vector<int>> newSubGraphs;
                            for (int &vertex: updatedSkylineCores[coreId].subGraphs[j]) {
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
                    //// reconstruct the subgraphs
                    vector<int> cores;
                    if (!start) {
                        bound = MultiBoundsQuery(nowPoint);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
                        start = true;
                    } else {
                        bound = MultiBoundsConnect(nowPoint);
                        originTree.QueryCore(RTree::AcceptEnclosing(bound), Visitor(), cores);
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

                    //// Calculate connectivity
                    for (int j = 0; j < tempResults.size(); j++) {
                        for (int node: updatedGraph.neighbors[tempResults[j]]) {
                            if (nodes[node]) {
                                UF.Unite(tempResults[j], node);
                            }
                        }
                    }
                    if (newCnt >= updatedTraverseVertexNum[i]) {
                        break;
                    }
                }
//                // updated R-tree node
                if (changedNodeCo.find(nowPoint) != changedNodeCo.end()) {
                    auto RNodes = updatedMapCoNode.find(nowPoint)->second;
                    vector<bool> tempNodes;
                    tempNodes.resize(vertexNum);
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
                    changedNodeCo.erase(nowPoint);
                }
            }
        }
        end_time = clock();
        total_time += (double) (end_time - start_time) / CLOCKS_PER_SEC;
    }
    cout << "Update time : " << total_time / (double) testCase.size() << "s" << endl;
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
        if (i != 0 && coreMax[traverseVec[i]] == coreMax[traverseVec[i - 1]]) {
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
        for (auto &element: preResult) {
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
                for (int node: multilayer_graph.neighbors[tempResults[j]]) {
                    if (nodes[node] && UF2.Find(tempResults[j]) != UF2.Find(node)) {
                        UF2.Unite(tempResults[j], node);
                        int tempCoreId = newVertexCore[node];
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

    dataset = "higgs";
//    dataset = "homo";
//  dataset = "sacchcere";
//    dataset = "dblp";
//  dataset = "example";

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

    // calculate index size of R*tree
    tree.indexSizeMulti(RTree::AcceptAny(), Visitor(), indexSize, mapCoreId);
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

    cout << "Index Size of R-tree and Skyline Cores: " << indexSize / 1024.0 / 1024.0 << "MB" << endl;
    int subgraphNum = 0;
    int connectEdgeNum = 0;
    for (auto &core: skylineCores) {
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

    cout << "Index Size of R-tree and Skyline Cores and Tags and Edges: " << indexSize / 1024.0 / 1024.0 << "MB"
         << endl;

    vector<int> nums = {1, 4, 16, 64, 256, 1024};
//    vector<int> nums = {1, 4};
    for (int i = 0; i < nums.size(); i++) {
        addOne(tree, multilayer_graph, nums[i]);
    }
    for (int i = 0; i < nums.size(); i++) {
        deleteOne(tree, multilayer_graph, nums[i]);
    }
    return 0;
}
