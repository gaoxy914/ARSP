#ifndef POINT_H
#define POINT_H

#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <deque>
#include <random>
#include <fstream>
#include <map>
#include <unordered_map>
#include <cstring>
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <queue>
#include <iostream>
extern "C" {
    #include "glpk.h"
}

// #define _DEBUG_

#define PREC 6
#define PI 3.14159265

using namespace std;

class Point {
public:
    int m_dim;
    double *m_coord;

    Point();
    Point(const int& dim, const double* coord = nullptr);
    Point(const Point& other);
    virtual ~Point();
    Point& operator =(const Point& other);
    double operator [](const int& index) const;
    
    friend bool operator ==(const Point& t, const Point& s);
    friend ostream & operator <<(ostream& out, const Point& p);
};

#endif