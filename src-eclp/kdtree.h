#ifndef KDTREE_H
#define KDTREE_H

#include "instance.h"
#include "hyperrect.h"

class KDTree {
    struct Node {
        int m_axis;
        int m_depth;
        int *m_index;
        int m_cnt;
        Node *m_left;
        Node *m_right;
        HyperRect m_mbr;
    };

    struct Iterator {
        int *m_sigma;
        double m_beta;
        int m_xi;
        vector<int> m_index;
    };

    int _dim;
    int _n;
    Instance *_data;
    Node *_root;
    Iterator _iter;
    int _m;

    void _clear(Node *node);
    void _init_iterator(const int& m);
    void _mix_sky_prob(int *index, int cnt, int depth, const int* obj_cnt, unordered_map<int, double>& result);

public:
    KDTree(const int& dim);
    virtual ~KDTree();

    void load_data(const int& n, const Instance* data, const vector<Point>& vertices = vector<Point>());
    void mix_sky_prob(const int& m, const int* obj_cnt, unordered_map<int, double>& result);
};

#endif