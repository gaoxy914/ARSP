#ifndef __OBJECT__
#define __OBJECT__

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <math.h>
#include <deque>
#include <random>
#include <fstream>
#include <map>
#include <unordered_map>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

#define _LINUX_

#ifdef _LINUX_
    #include <sys/time.h>
#endif

extern "C" {
    #include "glpk.h"
}

using namespace std;

inline double dot(const int& dim, const double *u, const double *v) {
    double sum = 0.0;
    for (int i = 0; i < dim; ++ i) sum += u[i]*v[i];
    return sum;
}

class HyperBox {
public:
    int dim;
    double *left_bottom, *right_top;

    HyperBox() {
        left_bottom = right_top = nullptr;
    }

    HyperBox(const int& dim) : dim(dim) {
        left_bottom = new double[dim];
        right_top = new double[dim];
    }

    HyperBox(const int& dim, double* coord) {
        this->dim = dim;
        this->left_bottom = new double[dim];
        this->right_top = new double[dim];
        for (int i = 0; i < dim; ++ i) {
            this->left_bottom[i] = coord[i];
            this->right_top[i] = coord[i];
        }
    }

    HyperBox(const int& dim, double* left_bottom, double* right_top) {
        this->dim = dim;
        this->left_bottom = new double[dim];
        this->right_top = new double[dim];
        for (int i = 0; i < dim; ++ i) {
            this->left_bottom[i] = left_bottom[i];
            this->right_top[i] = right_top[i];
        }
    }

    HyperBox(const HyperBox& box) {
        dim = box.dim;
        left_bottom = new double[dim];
        right_top = new double[dim];
        for (int i = 0; i < dim; ++ i) {
            left_bottom[i] = box.left_bottom[i];
            right_top[i] = box.right_top[i];
        }
    }

    ~HyperBox() {
        if (left_bottom) delete[] left_bottom;
        if (right_top) delete[] right_top;
    }

    /* return centroid of the box */
    double* GetCentroid() const {
        double *centroid = new double[dim];
        for (int i = 0; i < dim; ++ i) {
            centroid[i] = (left_bottom[i] + right_top[i])/2;
        }
        return centroid;
    }

    /* MBR of two box */
    friend HyperBox operator +(const HyperBox& _boxa, const HyperBox& _boxb) {
        HyperBox box(_boxa.dim, _boxa.left_bottom, _boxa.right_top);
        for (int i = 0; i < _boxa.dim; ++ i) {
            box.left_bottom[i] = min(box.left_bottom[i], _boxb.left_bottom[i]);
            box.right_top[i] = max(box.right_top[i], _boxb.right_top[i]);
        }
        return box;
    }

    void Merge(const double* coord) {
        for (int i = 0; i < dim; ++ i) {
            left_bottom[i] = min(coord[i], left_bottom[i]);
            right_top[i] = max(coord[i], right_top[i]);
        }
    }

    /* true : left_bottom == right_top */
    bool IsPoint() const {
        for (int i = 0; i < dim; ++ i) if (left_bottom[i] != right_top[i]) return false;
        return true;
    }

    /* true : right_top R-dominates instance */
    bool RDominates(const double* coord, const HyperBox& box) const {
        double sum = coord[dim - 1] - right_top[dim - 1];
        for (int i = 0; i < dim - 1; ++ i) {
            if (coord[i] - right_top[i] > 0) sum += (coord[i] - right_top[i])*box.left_bottom[i];
            else sum += (coord[i] - right_top[i])*box.right_top[i];
        }
        return sum >= 0;
    }

    HyperBox GetSubSpace(const int& k) const {
        HyperBox box(dim);
        for (int i = 0; i < dim; ++ i) {
            box.left_bottom[i] = left_bottom[i] + (right_top[i] - left_bottom[i])/2*((k>>i)&1);
            box.right_top[i] = box.left_bottom[i] + (right_top[i] - left_bottom[i])/2;
        }
        return box;
    }

    int PointLocation(const double* q_point, HyperBox& subspace) const {
        int k = 0;
        for (int i = 0; i < dim; ++ i) {
            if (q_point[i] > left_bottom[i] + (right_top[i] - left_bottom[i])/2) {
                k += (1<<i);
                subspace.left_bottom[i] = left_bottom[i] + (right_top[i] - left_bottom[i])/2;
            } else subspace.left_bottom[i] = left_bottom[i];
            subspace.right_top[i] = subspace.left_bottom[i] + (right_top[i] - left_bottom[i])/2;
        }
        return k;
    }

    friend ostream & operator <<(ostream& out, const HyperBox& box) {
        if (box.IsPoint()) {
            cout << "(";
            for (int i = 0; i < box.dim - 1; ++ i) cout << box.left_bottom[i] << ", ";
            cout << box.left_bottom[box.dim - 1] << ")";
            return out;
        }
        cout << "(";
        for (int i = 0; i < box.dim - 1; ++ i) cout << box.left_bottom[i] << ", ";
        cout << box.left_bottom[box.dim - 1] << ")\t";
        cout << "(";
        for (int i = 0; i < box.dim - 1; ++ i) cout << box.right_top[i] << ", ";
        cout << box.right_top[box.dim - 1] << ")";
        return out;
    }
};

/* format : x[d] = w[1]x[1] + ... + w[d-1]x[d-1] + w[d] */
class HyperPlane {
public:
    int dim;
    int obj_id;
    int ins_id;
    double prob;
    double *coef; // w[1], ..., w[d], i in [0, d - 2] : w[i] = coord[i], w[d - 1] = -coord[d - 1]

    HyperPlane() : obj_id(-1), ins_id(-1), prob(0), coef(nullptr) {}

    HyperPlane(const int& dim, const int& obj_id, const int& ins_id, const double& prob, const double* coef) {
        this->obj_id = obj_id;
        this->ins_id = ins_id;
        this->prob = prob;
        this->coef = new double[dim];
        for (int i = 0; i < dim; ++ i) this->coef[i] = coef[i];
    }

    HyperPlane(const HyperPlane& plane) {
        dim = plane.dim;
        obj_id = plane.obj_id;
        ins_id = plane.ins_id;
        prob = plane.prob;
        coef = new double[dim];
        for (int i = 0; i < dim; ++ i) coef[i] = plane.coef[i];
    }

    ~HyperPlane() {
        if (coef) delete[] coef;
    }

    /*
     * object funciton : z = x[d] - w[1]x[1] - ... - w[d-1]x[d-1]
     * subject to : x[i] \in [region.l[i], region.r[i]] for i \in {1, ..., d}
     */
    void CalExtremes(const HyperBox& region, double& max_value, double& min_value) const {
        max_value = region.right_top[dim - 1];
        min_value = region.left_bottom[dim - 1];
        for (int i = 0; i < dim - 1; ++ i) {
            max_value -= coef[i] > 0 ? coef[i]*region.left_bottom[i] : coef[i]*region.right_top[i];
            min_value -= coef[i] > 0 ? coef[i]*region.right_top[i] : coef[i]*region.left_bottom[i];
        }
    }

    bool Intersect(const HyperBox& region) const {
        double max_value = 0, min_value = 0;
        CalExtremes(region, max_value, min_value);
        return min_value < coef[dim - 1] && max_value > coef[dim - 1];
    }

    bool Above(const HyperBox& region) const {
        double max_value = 0, min_value = 0;
        CalExtremes(region, max_value, min_value);
        return max_value <= coef[dim - 1];
    }

    bool RDominates(const HyperPlane& plane, const HyperBox& R) const {
        double sum = -plane.coef[dim - 1] + coef[dim - 1];
        for (int i = 0; i < dim; ++ i ) {
            if (plane.coef[i] > coef[i]) sum += (plane.coef[i] - coef[i])*R.left_bottom[i];
            else sum += (plane.coef[i] - coef[i])*R.right_top[i];
        }
        return sum >= 0;
    }

    void GetQueryPoints(const HyperBox& R, vector<double*>& points) const {
        vector<double*> vertices;
        // R.GetVertices(vertices);
        points.reserve(vertices.size());
        for (int i = 0; i < vertices.size(); ++ i) {
            points[i] = new double[dim];
            points[i][dim - 1] = coef[dim - 1];
            for (int j = 0; j < dim - 1; ++ i) {
                points[i][j] = -vertices[i][j];
                points[i][dim - 1] -= vertices[i][j]*coef[j];
            }
        }
    }
};

class InstanceBase {
public:
    int dim;
    int obj_id;
    int ins_id;
    double prob; // precision concern, now for cnt in obj
    double *coord;
    
    InstanceBase() : obj_id(-1), ins_id(-1), prob(0), coord(nullptr) {}

    InstanceBase(const int& dim, const int& obj_id, const int& ins_id, const double& prob, const double* coord) {
        this->dim = dim;
        this->obj_id = obj_id;
        this->ins_id = ins_id;
        this->prob = prob;
        this->coord = new double[dim];
        for (int i = 0; i < dim; ++i) this->coord[i] = coord[i];
    }

    InstanceBase(const InstanceBase& instance) {
        dim = instance.dim;
        obj_id = instance.obj_id;
        ins_id = instance.ins_id;
        prob = instance.prob;
        coord = new double[dim];
        for (int i = 0; i < dim; ++ i) coord[i] = instance.coord[i];
    }

    InstanceBase& operator =(const InstanceBase& instance) {
        if (&instance != this) {
            dim = instance.dim;
            obj_id = instance.obj_id;
            ins_id = instance.ins_id;
            prob = instance.prob;
            for (int i = 0; i < dim; ++ i) coord[i] = instance.coord[i];
        }
        return *this;
    }

    ~InstanceBase() {
        if (coord) delete[] coord;
    }

    double Score(const vector<double>& weight) const {
        double sum = 0.0;
        for (int i = 0; i < dim; ++ i) sum += weight[i]*coord[i];
        return sum;
    }

    /* dominate test for loop based algorithm */
    bool Dominates(const InstanceBase& instance) const {
        if (ins_id == instance.ins_id) return false;
        for (int i = 0; i < dim; ++ i) {
            if (coord[i] > instance.coord[i]) return false;
        }
        return true;
    }

    /* R-dominate test */
    bool RDominates(const InstanceBase& instance, const HyperBox& box) const {
        double sum = instance.coord[dim - 1] - coord[dim - 1];
        for (int i = 0; i < dim - 1; ++ i) {
            if (instance.coord[i] - coord[i] > 0) sum += (instance.coord[i] - coord[i])*box.left_bottom[i];
            else sum += (instance.coord[i] - coord[i])*box.right_top[i];
        }
        return sum >= 0;
    }

    /* dominate test for transform based algorithm */
    bool Dominates(const double *point) const {
        bool equal = true;
        for (int i = 0; i < dim; ++ i) {
            if (coord[i] > point[i]) return false;
            else if (coord[i] < point[i]) equal = false;
        }
        return !equal;
    }
};

#endif