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

#define Dim 3 // data dimension

extern "C" {
    #include "glpk.h"
}

using namespace std;

inline double dot(double *u, double *v) {
    double sum = 0.0;
    for (int i = 0; i < Dim; ++ i) sum += u[i]*v[i];
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

    /* return 2^d vertices */
    void GetVertices(vector<double*>& vertices) const {
        vertices.resize(int(pow(2, dim)));
        for (int i = 0; i < pow(2, dim); ++ i) {
            int k = i;
            vertices[i] = new double[dim + 1];
            for (int j = 0; j < dim; ++ j) {
                vertices[i][j] = ((k>>j)&1) == 0 ? left_bottom[j] : right_top[j];
            }
            vertices[i][dim] = 1;
        }
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
    int obj_id;
    int ins_id;
    double prob;
    double *coef; // w[1], ..., w[d], i in [0, d - 2] : w[i] = coord[i], w[d - 1] = -coord[d - 1]

    HyperPlane() : obj_id(-1), ins_id(-1), prob(0), coef(nullptr) {}

    HyperPlane(const int& obj_id, const int& ins_id, const double& prob, const double* coef) {
        this->obj_id = obj_id;
        this->ins_id = ins_id;
        this->prob = prob;
        this->coef = new double[Dim];
        for (int i = 0; i < Dim; ++ i) this->coef[i] = coef[i];
    }

    HyperPlane(const HyperPlane& plane) {
        obj_id = plane.obj_id;
        ins_id = plane.ins_id;
        prob = plane.prob;
        coef = new double[Dim];
        for (int i = 0; i < Dim; ++ i) coef[i] = plane.coef[i];
    }

    /*
     * object funciton : z = x[d] - w[1]x[1] - ... - w[d-1]x[d-1]
     * subject to : x[i] \in [region.l[i], region.r[i]] for i \in {1, ..., d}
     */
    void CalExtremes(const HyperBox& region, double& max_value, double& min_value) const {
        max_value = region.right_top[Dim - 1];
        min_value = region.left_bottom[Dim - 1];
        for (int i = 0; i < Dim - 1; ++ i) {
            max_value -= coef[i] > 0 ? coef[i]*region.left_bottom[i] : coef[i]*region.right_top[i];
            min_value -= coef[i] > 0 ? coef[i]*region.right_top[i] : coef[i]*region.left_bottom[i];
        }
    }

    bool Intersect(const HyperBox& region) const {
        double max_value = 0, min_value = 0;
        CalExtremes(region, max_value, min_value);
        return min_value < coef[Dim - 1] && max_value > coef[Dim - 1];
    }

    bool Above(const HyperBox& region) const {
        double max_value = 0, min_value = 0;
        CalExtremes(region, max_value, min_value);
        return max_value <= coef[Dim - 1];
    }

    bool RDominates(const HyperPlane& plane, const HyperBox& R) const {
        double sum = -plane.coef[Dim - 1] + coef[Dim - 1];
        for (int i = 0; i < Dim; ++ i ) {
            if (plane.coef[i] > coef[i]) sum += (plane.coef[i] - coef[i])*R.left_bottom[i];
            else sum += (plane.coef[i] - coef[i])*R.right_top[i];
        }
        return sum >= 0;
    }

    void GetQueryPoints(const HyperBox& R, vector<double*>& points) const {
        vector<double*> vertices;
        R.GetVertices(vertices);
        points.reserve(vertices.size());
        for (int i = 0; i < vertices.size(); ++ i) {
            points[i] = new double[Dim];
            points[i][Dim - 1] = coef[Dim - 1];
            for (int j = 0; j < Dim - 1; ++ i) {
                points[i][j] = -vertices[i][j];
                points[i][Dim - 1] -= vertices[i][j]*coef[j];
            }
        }
    }
};

class InstanceBase {
public:
    int obj_id;
    int ins_id;
    double prob; // precision concern, now for cnt in obj
    double *coord;
    
    InstanceBase() : obj_id(-1), ins_id(-1), prob(0), coord(nullptr) {}

    InstanceBase(const int& obj_id, const int& ins_id, const double& prob, const double* coord) {
        this->obj_id = obj_id;
        this->ins_id = ins_id;
        this->prob = prob;
        this->coord = new double[Dim];
        for (int i = 0; i < Dim; ++i) this->coord[i] = coord[i];
    }

    InstanceBase(const InstanceBase& instance) {
        obj_id = instance.obj_id;
        ins_id = instance.ins_id;
        prob = instance.prob;
        coord = new double[Dim];
        for (int i = 0; i < Dim; ++ i) coord[i] = instance.coord[i];
    }

    /* dominate test for loop based algorithm */
    bool Dominates(const InstanceBase& instance) const {
        if (ins_id == instance.ins_id) return false;
        for (int i = 0; i < Dim; ++ i) {
            if (coord[i] > instance.coord[i]) return false;
        }
        return true;
    }

    /* R-dominate test */
    bool RDominates(const InstanceBase& instance, const HyperBox& box) const {
        double sum = instance.coord[Dim - 1] - coord[Dim - 1];
        for (int i = 0; i < Dim - 1; ++ i) {
            if (instance.coord[i] - coord[i] > 0) sum += (instance.coord[i] - coord[i])*box.left_bottom[i];
            else sum += (instance.coord[i] - coord[i])*box.right_top[i];
        }
        return sum >= 0;
    }

    /* dominate test for transform based algorithm */
    bool Dominates(const int& dim, const double *other) {
        int equal = true;
        for (int i = 0; i < dim; ++ i) {
            if (coord[i] > other[i]) return false;
            else if (coord[i] < other[i]) equal = false;
        }
        return !equal;
    }
};

#endif