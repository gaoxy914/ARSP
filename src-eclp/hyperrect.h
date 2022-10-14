#ifndef HYPERRECT_H
#define HYPERRECT_H

#include "point.h"

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
    void append(const vector<double>& point);
    void get_vertices(vector<Point>& vertices);
    bool upper_edominate(const HyperRect& range) const;
    bool lower_edominate(const HyperRect& range) const;

    friend HyperRect operator +(const HyperRect& rect1, const HyperRect& rect2);
    friend ostream & operator <<(ostream& out, const HyperRect& rect);
};

#endif