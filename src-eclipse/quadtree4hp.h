#ifndef QUADTREE_FOR_HYPERPLANE_H
#define QUADTREE_FOR_HYPERPLANE_H

#include "instance.h"
#include "hyperrect.h"

#define MAXPLANES 8
// #define MAXDEPTH 6

struct HyperPlane {
    Instance *m_pt;
    Instance *m_ps;
    vector<double> m_plane;
};

class QuadTree4HP {
public:
    struct Node {
        vector<Node*> _children;
        vector<int> _indices;
    };

    int _dim;
    int _n;
    int _nchild;
    vector<HyperPlane> _hyperplanes;
    Node *_root;
    HyperRect _space;

    Node *_build_rec(HyperRect& space, vector<int> indices, int depth);
    void _query_rec(Node *node, const HyperRect& rect, HyperRect& space, vector<int>& result);
    QuadTree4HP(const int& dim);
    virtual ~QuadTree4HP();
    void build(const vector<HyperPlane>& hyperplanes);
    void rect_query(const HyperRect& rect, vector<HyperPlane>& result);
};

#endif