#ifndef KDTREE_FOR_SHIFT_POINT_H
#define KDTREE_FOR_SHITT_POINT_H

#include "instance.h"
#include "hyperrect.h"

struct ShiftPoint {
    Instance *m_pt;
    Instance *m_ps;
    // vector<double> m_coord;

    double operator [](const int& index) const {
        return (m_pt->m_point[index] - m_ps->m_point[index]);
    }
};

class KDTree4SP {
public:
    struct Node {
        // int *m_index;
        // int m_cnt;
        // unordered_map<int, bool> m_domiante;
        vector<pair<int, int> > m_sigma;
        Node *m_left;
        Node *m_right;
        // HyperRect m_mbr;
        double *m_lower;
        double *m_upper;
    };

    int _dim;
    int _n;
    ShiftPoint *_data;
    Node *_root;
    Instance _origin;

    int _obj_cnt;

    int _ori_n;
    int *_sigma;

    bool _dominate(double *point, const HyperRect& rect);

    Node *_build_rec(int *index, int cnt, int depth);
    void _range_query_rec(const Node* node, const HyperRect& rect, unordered_map<int, double>& result);
    KDTree4SP(const int& dim, const int& ori_n, const int& obj_cnt);
    virtual ~KDTree4SP();
    void build(const int& n, ShiftPoint *data);
    void range_query(const HyperRect& rect, unordered_map<int, double>& result);
};

#endif