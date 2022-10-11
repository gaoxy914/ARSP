#ifndef ECLIPSE_H
#define ECLIPSE_H

#include "instance.h"
#include "quadtree4hp.h"
#include "kdtree4sp.h"

class Dataset {
public:
    int _dim;
    int _n;
    int _sky;
    Instance *_data;
    QuadTree4HP *_qd4hp;
    KDTree4SP *_kd4sp;

    Dataset(const int& dim, const int& n);
    virtual ~Dataset();

    void gen_data();
    void load_data();

    void build_quadtree();
    void build_multi_tree();

    void baseline(const HyperRect& rect, vector<int>& result);
    void quad(const HyperRect& rect, vector<int>& result);
    void multi(const HyperRect& rect, unordered_map<int, bool>& result);
};

#endif