#ifndef HYPERRECT_H
#define HYPERRECT_H

#include "point.h"

// for fast spherical volume
const double UNIT_SPHERE_VOLUMES[] = { // Dimension
    0.000000, 2.000000, 3.141593, //   0,1,2
    4.188790, 4.934802, 5.263789, //   3,4,5
    5.167713, 4.724766, 4.058712, //   6,7,8
    3.298509, 2.550164, 1.884104, //   9,10,11
    1.335263, 0.910629, 0.599265, //   12,13,14
    0.381443, 0.235331, 0.140981, //   15,16,17
    0.082146, 0.046622, 0.025807  //   18,19,20 
};

class HyperRect {
public:
    int m_dim;
    double *m_lower;
    double *m_upper;

    HyperRect();
    HyperRect(const int& dim, const double* lower = nullptr, const double* upper = nullptr);
    HyperRect(const Point& lower, const Point& upper);
    HyperRect(const HyperRect& other);
    virtual ~HyperRect();
    HyperRect& operator =(const HyperRect& other);

    void append(const Point& point);
    void append(const HyperRect& rect);
    double spherical_volume() const;
    double lower_score(const Point& weight) const;
    bool upper_fdominate_V(const Point& point, const vector<Point>& vertex) const;
    bool lower_fdominate_V(const Point& point, const vector<Point>& vertex) const;
    bool upper_dominate(const Point& point) const;
    bool lower_dominate(const Point& point) const;
    void z_split(vector<HyperRect>& rect, Point& mid_point);
    HyperRect* get_subrect(const int& k);
    Point get_mid();

    friend HyperRect operator +(const HyperRect& rect1, const HyperRect& rect2);
    friend ostream & operator <<(ostream& out, const HyperRect& rect);
};

#endif