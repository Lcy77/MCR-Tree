#ifndef RSTARTREE_UNIONFIND_H
#define RSTARTREE_UNIONFIND_H

#include <vector>
#include <unordered_set>
#include <set>
#include <map>

using namespace std;

class UnionFind {
public:
    int node_num;
    vector<int> father;
    vector<int> rank;
    UnionFind(int n);
    void Init();
    void Init(vector<int> &nums);
    int Find(int x);
    void Unite(int x, int y);
    void UniteRecord(int x, int y, vector<int> &revisedVertex);
    int UniteRoot(int x, int y);
    void BatchUnite(set<int> &nodes);
};
#endif //RSTARTREE_UNIONFIND_H
