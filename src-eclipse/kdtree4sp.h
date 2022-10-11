#ifndef KDTREE_FOR_SHIFT_POINT_H
#define KDTREE_FOR_SHITT_POINT_H

#include "instance.h"
#include "hyperrect.h"

#define MAXPLANES 8

struct ShiftPoint {
    Instance *m_pt;
    Instance *m_ps;
    vector<double> m_coord;
};

class KDTree4SP {
public:
    struct Node {
        int *m_index;
        int m_cnt;
        unordered_map<int, bool> m_domiante;
        Node *m_left;
        Node *m_right;
        HyperRect m_mbr;
    };

    int _dim;
    int _n;
    ShiftPoint *_data;
    Node *_root;
    Instance _origin;

    Node *_build_rec(int *index, int cnt, int depth);
    void _range_query_rec(const Node* node, const HyperRect& rect, unordered_map<int, bool>& result);
    KDTree4SP(const int& dim);
    virtual ~KDTree4SP();
    void build(const int& n, ShiftPoint *data);
    void range_query(const HyperRect& rect, unordered_map<int, bool>& result);
};

#endif