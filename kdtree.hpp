#ifndef __KDTREE__
#define __KDTREE__

#include "object.hpp"

class KDTree {
    struct Node {
        int axis; // split dimension's axis
        int* indices;
        int cnt;
        Node *left, *right;
        HyperBox R; // mbr of this node
        Node() : axis(-1), cnt(0), indices(nullptr), left(nullptr), right(nullptr) {}
    };

    Node *root; // root node;
    int dim;
    vector<InstanceBase> points; // points

    Node* BuildRecursive(int* indices, int npoints, int depth) {
        if (npoints <= 0) return nullptr;
        const int axis = depth%dim;
        const int mid = (npoints - 1)/2;
        nth_element(indices, indices + mid, indices + npoints, [&](int lhs, int rhs){
            return points[lhs][axis] < points[rhs][axis];
        });
        Node *node = new Node();
        node->axis = axis;
        node->cnt = npoints;
        node->indices = new int[npoints];
        for (int i = 0; i < npoints; ++ i) {
            node->indices[i] = indices[i];
            // node->R = node->R + points[indices[i]].coord;
        }
        node->left = BuildRecursive(indices, mid, depth + 1);
        node->right = BuildRecursive(indices + mid, npoints - mid, depth + 1);
    }

    void GetProbRecursive(map<int, double> results) const {

    }
public:
};

#endif