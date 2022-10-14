#ifndef INSTANCE_H
#define INSTANCE_H

#include "hyperrect.h"

class Instance {
public:
    Point m_point;
    int m_id;
    int m_obj_id;

    Instance();
    Instance(const int& dim, const double* coord = nullptr, const int& id = -1, const int& obj_id = -1);
    Instance(const Instance& other);
    virtual ~Instance();
    double score(const Point& weight) const;
    bool dominate(const HyperRect& rect) const;
    bool dominate(const Instance& other) const;
    bool cross(const HyperRect& rect) const;
    Instance& operator =(const Instance& other);

    friend ostream & operator <<(ostream& out, const Instance& t);

    struct Comparator {
        Point m_weight;
        Comparator(const Point& weight) : m_weight(weight) {}

        bool operator ()(const Instance& t, const Instance& s) const {
            return t.score(m_weight) < s.score(m_weight);
        }
    };

};

#endif