#ifndef QUADTREE_H
#define QUADTREE_H

#include "instance.h"
#include "hyperrect.h"
#include "region.h"

class QuadTree {
    struct Node {
        int m_order;
        // vector<int> m_index;
        // Node **m_children;
        HyperRect m_mbr;
    };

    struct Iterator {
        int *m_sigma;
        double m_beta;
        int m_xi;
        vector<int> m_index;
        vector<int> tmp;
    };

    int _dim;
    int _n;
    int _nchild;
    vector<Instance> _data;
    Node *_root;
    Iterator _iter;

    void _init_iterator(const int& m);
    void _clear(Node *node);
    void _mix_sky_prob(vector<int> index, const int& depth, const int* obj_cnt, unordered_map<int, double>& result);

    // Node *_build(vector<int> index, const int& depth);
    // void _sky_porb(const Node* node, const int* obj_cnt, unordered_map<int, double>& result);

public:
    QuadTree(const int& dim);
    virtual ~QuadTree();

    void load_data(const int& n, const Instance* data, Region* weights = nullptr);
    // void build();
    // void sky_prob(const int& m, const int* obj_cnt, unordered_map<int, double>& result);
    void mix_sky_prob(const int& m, const int* obj_cnt, unordered_map<int, double>& result);
};

#endif