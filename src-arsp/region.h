#ifndef REGION_H
#define REGION_H

#include "point.h"

#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/QhullUser.h"
#include "libqhullcpp/QhullVertex.h"

using orgQhull::Coordinates;
using orgQhull::Qhull;
using orgQhull::QhullFacetListIterator;
using orgQhull::QhullHyperplane;

class HalfSpace {
public:
    int m_dim; // number of variables
    double* m_coef; // a, b (length = m_dim + 1)
    bool m_side; // true: a x w >= b

    HalfSpace();
    HalfSpace(const int& dim, const double* coef = nullptr, const bool& side = false);
    HalfSpace(const HalfSpace& other);
    virtual ~HalfSpace();
    HalfSpace& operator =(const HalfSpace& other);

    friend ostream & operator <<(ostream& out, const HalfSpace& h);
};

/*
 * Region: linear constraints A x w >= b
 * intersection of halfspaces in R^{d-1} since w[d - 1] is replaced by 1 - w[0] ... - w[d - 2]
 * m_boundplane R^{d-1} : bounding hyperplane of halfspace
 * m_vertex in R^d : vertices of region,
 * m_inner in R^d: inner point of region
 */
class Region {

    int _dim;
    int _c;
    vector<HalfSpace> _boundplane;
    vector<Point> _vertex;
    Point _inner;
    // string _qpath_prefix;

    glp_prob *_lp;

    double _roundoff(const double& value, unsigned char prec);
    double _rand_uniform(const double& a, const double& b);
    double _rand_normal(const double& med, const double& var);

public:
    Region(const int& dim, const int& c_cnt);
    virtual ~Region();

    /* bool compute_inner(const int& dim, double *inner); */
    Point& get_inner();
    Point& get_vertex(const int& index);
    glp_prob* get_lp();
    vector<Point>& get_vertex();
    int get_vertex_size();
    /* 
    * generate linear constraints
    * 1. generate an inner point i in S^{d-1}
    * 2. pick a pivot p between some vertex of S^{d-1} and i
    * 3. pick a random vector a to form the hyperplane a(w - p) = 0
    * 4. choose the halfspace bounding by a(w - p) = 0 that contains i
    */
    void gen_query();
    void gen_weak_rankings();
    void load_query();
    void load_weak_rankings();
    void load_nba_query();
    void write_query();
    void add_plane();
    void build_lp();
    void compute_vertex();
    void print();
};

#endif