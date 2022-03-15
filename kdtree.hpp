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

    struct Info {
        vector<InstanceBase> cap_instances;
        double* sigma;
        double beta;
        int xi;
    };

    Node *root; // root node;
    int dim;
    vector<InstanceBase> points; // points

    Node* BuildRecursive(int* indices, int npoints, int depth) {
        if (npoints <= 0) return nullptr;
        const int axis = depth%dim;
        const int mid = (npoints - 1)/2;
        nth_element(indices, indices + mid, indices + npoints, [&](int lhs, int rhs){
            return points[lhs].coord[axis] < points[rhs].coord[axis];
        });
        Node *node = new Node();
        node->axis = axis;
        node->cnt = npoints;
        node->indices = new int[npoints];
        node->R = HyperBox(dim, points[*indices].coord);
        for (int i = 0; i < npoints; ++ i) {
            node->indices[i] = indices[i];
            node->R = node->R + HyperBox(dim, points[indices[i]].coord);
        }
        node->left = BuildRecursive(indices, mid, depth + 1);
        node->right = BuildRecursive(indices + mid, npoints - mid, depth + 1);
        return node;
    }

    void ClearRecursive(Node *node) {
        if (node == nullptr) return;
        if (node->left) ClearRecursive(node->left);
        if (node->right) ClearRecursive(node->right);
        delete node;
    }

    void CalSkyProbRecursive(const Node* node, Info& info, map<int, double>& results) const {
        // update info
        vector<InstanceBase> instances_pv = info.cap_instances;
        vector<InstanceBase> instances_v, domineats_v;
        for (auto ins : instances_pv) {
            if (ins.Dominates(dim, node->R.right_top)) { // dominates region
                domineats_v.push_back(ins);
            } else if (ins.Dominates(dim, node->R.left_bottom)) { // intersects with region
                instances_v.push_back(ins);
            }
        }
        if (node->cnt == 1) {
            int index = *(node->indices);
            int ins_id = points[index].ins_id, obj_id = points[index].obj_id;
            if (info.xi > 0) results[ins_id] = 0;
            else results[ins_id] = info.beta/(1 - info.sigma[obj_id]);
            // undo the update to info 
            for (auto ins : instances_v) {

            }
            for (auto ins : domineats_v) {
                
            }
        } else {
            if (node->left) {
                info.cap_instances = instances_v;
                CalSkyProbRecursive(node->left, info, results);
                // undo the update to info
                info.cap_instances = instances_pv;
                for (auto ins : instances_v) {

                }
                for (auto ins : domineats_v) {

                }
            }
            if (node->right) {
                info.cap_instances = instances_v;
                CalSkyProbRecursive(node->right, info, results);
                // undo the update to infor
                info.cap_instances = instances_pv;
                for (auto ins : instances_v) {

                }
                for (auto ins : domineats_v) {
                    
                }

            }
        }
    }

public:
    KDTree() : root(nullptr), dim(0) {};
    KDTree(const int& dim, const vector<InstanceBase>& points) : root(nullptr), dim(dim) { Build(points); }
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

    }
};

#endif