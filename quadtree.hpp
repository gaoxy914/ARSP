#ifndef __QUADTREE__
#define __QUADTREE__

#include "object.hpp"

class QuadTree {
    int m;
    int c;
    int n_child;

    struct Info {
        int *sigma;
        double beta;
        int xi;

        Info() : sigma(nullptr), beta(1), xi(0) {}
    };

    /*
     * internal node :
     * pointers to children
     * leaf node :
     * 1) planes through its region
     * 2) tree of planes above its region or info
     */
    struct Node {
        vector<Node*> children;
        int n_through;
        int* through_planes; // planes through the leaf
        void* tree_or_info;

        Node() : n_through(0), through_planes(nullptr), tree_or_info(nullptr) {}
    };

    Node *root;
    vector<HyperPlane> planes;
    HyperBox space;

    Node* BuildRecursive(const HyperBox& space, vector<int> above, vector<int> through, int depth, int level) {
        Node *node = new Node();
        if (through.size() <= c) {
            node->n_through = through.size();
            node->through_planes = new int[through.size()];
            for (int i = 0; i < through.size(); ++ i) node->through_planes[i] = through[i];
            if (level == 0) {
                // build information
                Info info;
                info.sigma = new int[m];
                memset(info.sigma, 0, m*sizeof(info.sigma));
                for (int p : above) {
                    if (info.sigma[planes[p].obj_id] + 1 == planes[p].prob) {
                        info.beta *= planes[p].prob;
                        info.xi += 1;
                        info.sigma[planes[p].obj_id] += 1;
                    } else {
                        info.beta *= (planes[p].prob - info.sigma[planes[p].obj_id] - 1)/(planes[p].prob - info.sigma[planes[p].obj_id]);
                        info.sigma[planes[p].obj_id] += 1;
                    }
                }
                node->tree_or_info = &info;
            } else {
                // build next level tree
                through.clear();
                node->tree_or_info = BuildRecursive(this->space, through, above, 1, level - 1);
            }
        } else {
            node->children.reserve(n_child);
            for (int i = 0; i < n_child; ++ i) {
                HyperBox subspace = space.GetSubSpace(i);
                vector<int> subthrough;
                vector<int> subabove;
                for (int p : through) {
                    if (planes[p].Above(subspace)) subabove.push_back(p);
                    else if (planes[p].Intersect(subspace)) subthrough.push_back(p);
                }
                node->children[i] = BuildRecursive(subspace, subabove, subthrough, depth + 1, level);
            }
        }
        return node;
    }

    void ClearRecursive(Node *node) {
        
    }

public:
    QuadTree() : m(0), root(nullptr), n_child(pow(2, Dim)) {}

    QuadTree(const int& m, const vector<HyperPlane>& planes, const HyperBox& space) {
        this->m = m;
        n_child = pow(2, Dim);
        root = nullptr;
        Build(planes, space);
    }

    ~QuadTree() {
        Clear();
    }

    void Clear() {
        ClearRecursive(root);
        root = nullptr;
    }

    void Build(const vector<HyperPlane>& planes, const HyperBox& space) {
        if (root != nullptr) Clear();

    }

    double CalProb(const HyperBox& R, const InstanceBase& instance) const {

    }

};

#endif