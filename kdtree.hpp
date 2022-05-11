#ifndef __KDTREE__
#define __KDTREE__

#include "object.hpp"

class KDTree {
    struct Node {
        int axis; // split dimension's axis
        int* indices;
        int cnt;
        Node *left, *right;
        HyperBox* R; // mbr of this node
        Node() : axis(-1), cnt(0), indices(nullptr), left(nullptr), right(nullptr) {}
    };

    struct Info {
        int *sigma;
        double beta;
        int xi;
        vector<int> cap_instances;
    };

    Node *root; // root node;
    int dim;
    int m;
    vector<InstanceBase> points; // points

    Node* BuildRecursive(int* indices, int npoints, int depth) {
        if (npoints <= 0) return nullptr;
        const int axis = depth%dim;
        const int mid = npoints/2;
        // cout << "depth: " << depth << " axis: " << axis << " points number: " << npoints << " mid: " << mid <<endl;
        nth_element(indices, indices + mid, indices + npoints, [&](int lhs, int rhs){
            return points[lhs].coord[axis] < points[rhs].coord[axis];
        });
        // cout << points[indices[mid - 1]].coord[axis] << ", " << points[indices[mid]].coord[axis] << ", " << points[indices[mid + 1]].coord[axis] << endl;
        Node *node = new Node();
        node->axis = axis;
        node->cnt = npoints;
        node->indices = new int[npoints];
        node->R = new HyperBox(dim, points[*indices].coord);
        for (int i = 0; i < npoints; ++ i) {
            node->indices[i] = indices[i];
            // *node->R = *node->R + HyperBox(dim, points[indices[i]].coord);
            node->R->Merge(points[indices[i]].coord);
        }
        // cout << *node->R << endl;
        if (npoints > 1) {
            node->left = BuildRecursive(indices, mid, depth + 1);
            node->right = BuildRecursive(indices + mid, npoints - mid, depth + 1);
        }
        return node;
    }

    void ClearRecursive(Node *node) {
        if (node == nullptr) return;
        if (node->left) ClearRecursive(node->left);
        if (node->right) ClearRecursive(node->right);
        delete node;
    }

    void CalSkyProbRecursive(const Node* node, Info info, map<int, double>& results) const {
        // update info
        // cout << node->cnt << endl;
        vector<int> instances_pv = info.cap_instances;
        vector<int> instances_v, dominates_v;
        for (auto ins : instances_pv) {
            if (points[ins].Dominates(dim, node->R->left_bottom)) { // dominates region
                dominates_v.push_back(ins);
                if (info.sigma[points[ins].obj_id] + 1 == points[ins].prob) {
                    info.beta *= points[ins].prob;
                    info.xi += 1;
                    info.sigma[points[ins].obj_id] += 1;
                } else {
                    info.beta *= (points[ins].prob - info.sigma[points[ins].obj_id] - 1)/(points[ins].prob - info.sigma[points[ins].obj_id]);
                    info.sigma[points[ins].obj_id] += 1;
                }
            } else if (!node->R->IsPoint() && points[ins].Dominates(dim, node->R->right_top)) // intersects with region
                instances_v.push_back(ins);
        }
        info.cap_instances = instances_v;
        if (node->cnt == 1) {
            int index = *(node->indices);
            int ins_id = points[index].ins_id, obj_id = points[index].obj_id, prob = points[index].prob;
            if (info.xi > 0) results[ins_id] = 0;
            else results[ins_id] = (1/points[index].prob)*info.beta*prob/(prob - info.sigma[obj_id]);
            // cout << ins_id << '\t' << results[ins_id] << endl;
        } else {
            if (node->left) CalSkyProbRecursive(node->left, info, results);
            if (node->right) CalSkyProbRecursive(node->right, info, results);
        }
        for (auto ins : dominates_v) info.sigma[points[ins].obj_id] -= 1;
    }

public:
    KDTree() : root(nullptr), dim(0) {};
    KDTree(const int& dim, const int& m, const vector<InstanceBase>& points) : root(nullptr), dim(dim), m(m) { Build(points); }
    ~KDTree() { Clear(); }

    void Build(const vector<InstanceBase>& points) {
        if (root != nullptr) Clear();
        this->points = points;
        vector<int> indices(points.size());
        iota(begin(indices), end(indices), 0);
        root = BuildRecursive(indices.data(), (int)points.size(), 0);
    }

    void Clear() {
        ClearRecursive(root);
        root = nullptr;
        points.clear();
    }

    map<int, double> CalSkyPorb() const {
        map<int, double> results;
        Info info;
        info.sigma = new int[m];
        memset(info.sigma, 0, m*sizeof(int));
        // for (int i = 0; i < m; ++ i) info.sigma[i] = 0;
        info.beta = 1;
        info.xi = 0;
        // info.cap_instances = points;
        info.cap_instances.resize(points.size());
        iota(begin(info.cap_instances), end(info.cap_instances), 0);
        CalSkyProbRecursive(root, info, results);
        delete[] info.sigma;
        info.sigma = nullptr;
        return results;
    }
};

#endif