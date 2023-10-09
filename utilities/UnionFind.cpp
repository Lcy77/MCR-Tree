#include "UnionFind.h"

UnionFind::UnionFind(int n) {
    node_num = n;
    father.resize(node_num);
    rank.resize(node_num);
    for (int i = 0; i < node_num; i++) {
        father[i] = i;
        rank[i] = 0;
    }
}

void UnionFind::Init() {
    int i = 0;
    for (int& element : father) {
        element = i++;
    }
    for (int& element : rank) {
        element = 0;
    }
}

void UnionFind::Init(vector<int> &nums) {
    for (auto & element : nums) {
        father[element] = element;
        rank[element] = 0;
    }
}

int UnionFind::Find(int x) {
    if(x == father[x])
    return x;
    return father[x] = Find(father[x]); //在第一次查找时，将节点直连到根节点
}

void UnionFind::Unite(int x, int y){
    //查找到x，y的根节点
    x = Find(x);
    y = Find(y);
    if(x == y)
        return ;
    //判断两棵树的高度，然后在决定谁为子树
    if(rank[x] < rank[y]){
        father[x] = y;
    }else{
        father[y] = x;
        if(rank[x] == rank[y])
            rank[x]++;
    }
}

void UnionFind::UniteRecord(int x, int y, vector<int> &revisedVector) {
    //查找到x，y的根节点
    x = Find(x);
    y = Find(y);
    if(x == y)
        return ;
    //判断两棵树的高度，然后在决定谁为子树
    if(rank[x] < rank[y]){
        father[x] = y;
        revisedVector.push_back(x);
    }else{
        father[y] = x;
        revisedVector.push_back(y);
//        if(rank[x] == rank[y])
//            rank[x]++;
    }
}

int UnionFind::UniteRoot(int x, int y) {
    if (x == y)
        return x;
    if (rank[x] < rank[y]) {
        father[x] = y;
        return y;
    } else {
        father[y] = x;
        if(rank[x] == rank[y])
            rank[x]++;
        return x;
    }
}

//void UnionFind
void UnionFind::BatchUnite(set<int> &root_nodes) {
    if (root_nodes.size() == 0)
        return;
    vector<int> vec_root_nodes;
    for (int node : root_nodes) {
        vec_root_nodes.push_back(node);
    }
    for (int i = 0; i < vec_root_nodes.size() - 1; i++) {
        int x = vec_root_nodes[i];
        int y = vec_root_nodes[i+1];
        if(rank[x] < rank[y]){
            father[x] = y;
        }else{
            father[y] = x;
            if(rank[x] == rank[y])
                rank[x]++;
        }
    }
}
