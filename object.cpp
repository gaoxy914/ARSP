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

#define Dim 3

extern "C" {
    #include "glpk.h"
}

using namespace std;

inline double dot(const int& d, double *u, double *v) {
    double sum = 0.0;
    for (int i = 0; i < d; ++ i) sum += u[i]*v[i];
    return sum;
}

class HyperBox {
public:
    int d;
    double *left_bottom, *right_top;

    HyperBox() {
        left_bottom = right_top = nullptr;
    }

    HyperBox(const int& d, double* left_bottom, double* right_top) {
        this->d = d;
        this->left_bottom = new double[d];
        this->right_top = new double[d];
        for (int i = 0; i < d; ++ i) {
            this->left_bottom[i] = left_bottom[i];
            this->right_top[i] = right_top[i];
        }
    }
    
    ~HyperBox() {
        delete[] left_bottom;
        delete[] right_top;
        left_bottom = right_top = nullptr;
    }

    /* return centroid of the box */
    double* GetCentroid() const {
        double *centroid = new double[d];
        for (int i = 0; i < d; ++ i) {
            centroid[i] = (left_bottom[i] + right_top[i])/2;
        }
        return centroid;
    }

    /* return 2^d extremes */
    void GetExtremes(vector<double*> extremes) const {
        extremes.reserve(pow(2, d));
        for (int i = 1, k = 1; i < pow(2, d) - 1; ++ i, ++ k) {
            double *p = new double[d + 1];
            for (int j = 0; j < d; ++ j) {
                p[j] = ((k>>j)&1) == 0 ? left_bottom[j] : right_top[j];
            }
            p[d] = 1;
            extremes[i] = p;
        }
        extremes[0] = left_bottom;
        extremes[-1] = right_top;
    }

    /* Contain Test */
    bool Contain(const HyperBox& other) const {
        for (int i = 0; i < d; ++ i)
            if (left_bottom[i] > other.left_bottom[i] || right_top[i] < other.right_top[i]) return false;
        return true;
    }

    /* Intesection Test */
    bool Intersect(const HyperBox& other) const {
        for (int i = 0; i < d; ++ i)
            if (left_bottom[i] > other.right_top[i] || right_top[i] < other.left_bottom[i]) return false;
        
        return true;
    }

    /* return volumn */
    double GetArea() const {
        double area = 1;
        for (int i = 0; i < d; ++ i) area *= (right_top[i] - left_bottom[i]);
        return area;
    }

    /* return perimeter */
    double GetPerimeter() const {
        double perimeter = 0;
        for (int i = 0; i < d; ++ i) perimeter += (right_top[i] - left_bottom[i]);
        return perimeter;
    }

    /* return volumn of overlap */
    double GetOverlapArea(const HyperBox& other) const {
        double overlap = 1;
        if (!Intersect(other)) return 0;
        for (int i = 0; i < d; ++ i) {
            if (left_bottom[i] <= other.left_bottom[i]) {
                if (right_top[i] <= other.right_top[i]) overlap *= (right_top[i] - other.left_bottom[i]);
                else overlap *= (other.right_top[i] - other.left_bottom[i]);
            } else {
                if (right_top[i] <= other.right_top[i]) overlap *= (right_top[i] - left_bottom[i]);
                else overlap *= (other.right_top[i] - left_bottom[i]);
            }
        }
        return overlap;
    }

    /* MBR of two box */
    friend HyperBox operator +(const HyperBox& _boxa, const HyperBox& _boxb) {
        HyperBox box(_boxa.d, _boxa.left_bottom, _boxa.right_top);
        for (int i = 0; i < _boxa.d; ++ i) {
            box.left_bottom[i] = min(box.left_bottom[i], _boxb.left_bottom[i]);
            box.right_top[i] = max(box.right_top[i], _boxb.right_top[i]);
        }
        return box;
    }

    /* true : left_bottom == right_top */
    bool IsPoint() const {
        for (int i = 0; i < d; ++ i) if (left_bottom[i] != right_top[i]) return false;
        return true;
    }

    friend ostream & operator <<(ostream& out, const HyperBox& box) {
        if (box.IsPoint()) {
            cout << "(";
            for (int i = 0; i < box.d - 1; ++ i) cout << box.left_bottom[i] << ", ";
            cout << box.left_bottom[box.d - 1] << ")";
            return out;
        }
        cout << "(";
        for (int i = 0; i < box.d - 1; ++ i) cout << box.left_bottom[i] << ", ";
        cout << box.left_bottom[box.d - 1] << ")\t";
        cout << "(";
        for (int i = 0; i < box.d - 1; ++ i) cout << box.right_top[i] << ", ";
        cout << box.right_top[box.d - 1] << ")";
        return out;
    }
};

class Instance {
public:
    int obj_id;
    int ins_id;
    double prob;
    double *coord;
    double beta;
    // double *sigma;
    
    Instance() : obj_id(-1), ins_id(-1), prob(0), coord(nullptr), beta(1) {}

    Instance(const int& obj_id, const int& ins_id, const double& prob, const double* coord) {
        this->obj_id = obj_id;
        this->ins_id = ins_id;
        this->prob = prob;
        this->coord = new double[Dim];
        for (int i = 0; i < Dim; ++i) this->coord[i] = coord[i];
        this->beta = 1;
    }

    ~Instance() {
        if (coord != nullptr) delete[] coord;
        coord = nullptr;
    }

    /* return final result */
    double GetProb() const { return prob*beta; }

    /* dominate test */
    bool Dominates(const Instance& instance) const {
        for (int i = 0; i < Dim; ++ i) if (coord[i] > instance.coord[i]) return false;
        return true;
    }

    /* R-dominate test */
    bool RDominates(const Instance& instance, const HyperBox& box) const {
        double sum = instance.coord[Dim - 1] - coord[Dim - 1];
        for (int i = 0; i < Dim - 1; ++ i) {
            if (instance.coord[i] - coord[i] > 0) sum += (instance.coord[i] - coord[i])*box.left_bottom[i];
            else sum += (instance.coord[i] - coord[i])*box.right_top[i];
        }
        return sum >= 0;
    }
};

struct InstanceComparator {
    InstanceComparator(int _d, double *_weight) : d(_d), weight(_weight) {}
    bool operator ()(const Instance& u, const Instance& v) const {
        return dot(d, u.coord, weight) < dot(d, v.coord, weight);
    }
    double *weight;
    int d;
};

class DataOperator {
    int m, n, cnt;
    int* ins_cnt;
    vector<Instance*> dataset;
public:
    DataOperator() : m(0), n(0), cnt(0) {}

    DataOperator(const int& m, const int& cnt) {
        this->m = m;
        this->n = 0;
        this->cnt = cnt;
        this->ins_cnt = new int[m];
        dataset.resize(m);
    }

    ~DataOperator() {
        delete[] ins_cnt;
        for (auto iter : dataset) if (iter != nullptr) { delete[] iter; iter = nullptr; }
    }

    void PrintData() {
        cout << "m = " << m << "\t n = " << n << endl;
        for (int i = 0; i < m; ++ i) {
            cout << "obj " << i << "\t cnt = " << ins_cnt[i] << endl;
            for (int j = 0; j < ins_cnt[i]; ++ j) {
                cout << "(";
                for (int k = 0; k < Dim - 1; ++ k) cout << dataset[i][j].coord[k] << ", ";
                cout << dataset[i][j].coord[Dim - 1] << ")\t";
                cout << "prob = " << dataset[i][j].prob << endl;
            }
        }
    }

    void LoadData(const char *path) {
        ifstream file((string(path) + to_string(m) + string("/cnt.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) { file >> ins_cnt[i]; n += ins_cnt[i]; }
        file.close();
        file.open((string(path) + to_string(m) + string("/instances.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) {
            dataset[i] = new Instance[ins_cnt[i]];
            for (int j = 0; j < ins_cnt[i]; ++ j) {
                file >> dataset[i][j].obj_id;
                file >> dataset[i][j].ins_id;
                file >> dataset[i][j].prob;
                dataset[i][j].coord = new double[Dim];
                for (int k = 0; k < Dim; ++ k) file >> dataset[i][j].coord[k];
            }
        }
        file.close();
    }

    void WriteData(const char *path) {
        ofstream file((string(path) + to_string(m) + string("/cnt.data")).c_str(), ios::out);
        for (int i = 0; i < m; ++ i) file << ins_cnt[i] << " ";
        file.close();
        file.open((string(path) + to_string(m) + string("/instances.data")).c_str(), ios::out);
        for (int i = 0; i < m; ++ i) {
            for (int j = 0; j < ins_cnt[i]; ++ j) {
                file << dataset[i][j].obj_id << " ";
                file << dataset[i][j].ins_id << " ";
                file << dataset[i][j].prob << " ";
                for (int k = 0; k < Dim; ++ k) file << dataset[i][j].coord[k] << " ";
            }
        }
        file.close();
    }

    void GenData(const int& mode, const double& var1 = 0.05, const double& var2 = 0.025) {
        srand((unsigned)time(nullptr));
        switch (mode) {
        case 1:
            for (int i = 0; i < m; ++ i) {
                ins_cnt[i] = GenIndePoints(i, dataset[i]);
                n += ins_cnt[i];
            }
            break;
        case 2:
            for (int i = 0; i < m; ++ i) {
                ins_cnt[i] = GenAntiPoints(i, dataset[i], var1);
                n += ins_cnt[i];
            }
            break;
        case 3:
            for (int i = 0; i < m; ++ i) {
                ins_cnt[i] = GenCorrPoints(i, dataset[i], var1, var2);
                n += ins_cnt[i];
            }
            break;
        default:
            break;
        }
    }
private:
    double RandNormal(const double& med, const double& var) const {
        default_random_engine generator;
        normal_distribution<double> distribution(med, var);
        return distribution(generator);
    }

    /* generate independent instances with uniform distribution */
    int GenIndePoints(const int& obj_id, Instance* &points) const {
        double *center = new double[Dim];
        for (int i = 0; i < Dim; ++ i) center[i] = rand()/double(RAND_MAX);
        int ins_cnt = rand()%cnt + 1;
        double *l = new double[Dim];
        double *r = new double[Dim];
        for (int i = 0; i < Dim; ++ i) {
            double length = RandNormal(1.0, 0.025);
            l[i] = center[i] - length/2 > 0 ? center[i] - length/2 : 0;
            r[i] = center[i] + length/2 < 1 ? center[i] + length/2 : 1;
        }
        points = new Instance[ins_cnt];
        for (int i = 0; i < ins_cnt; ++ i) {
            points[i].obj_id = obj_id;
            points[i].ins_id = i;
            points[i].coord = new double[Dim];
            for (int j = 0; j < Dim; ++ j) points[i].coord[j] = rand()/double (RAND_MAX)*(r[j] - l[j]) + l[j];
            points[i].prob = 1/double(ins_cnt);
        }
        delete[] center;
        delete[] l;
        delete[] r;
        return ins_cnt;
    }

    /* generate anti-correlated instances with uniform distribution */
    int GenAntiPoints(const int& obj_id, Instance* &points, const double& var) const {
        double c = 0;
        do {c = RandNormal(0.5, var); } while (c >= 1 || c <= 0);
        double *center = new double[Dim];
        do {
            center[Dim - 1] = Dim*c;
            for (int i = 0; i < Dim - 1; ++ i) {
                center[i] = rand()/double(RAND_MAX);
                center[Dim - 1] -= center[i];
            }
        } while (center[Dim - 1] < 0 || center[Dim - 1] > 1);
        int ins_cnt = rand()%cnt + 1;
        double *l = new double[Dim];
        double *r = new double[Dim];
        for (int i = 0; i < Dim; ++ i) {
            double length = RandNormal(1.0, 0.025);
            l[i] = center[i] - length/2 > 0 ? center[i] - length/2 : 0;
            r[i] = center[i] + length/2 < 1 ? center[i] + length/2 : 1;
        }
        points = new Instance[ins_cnt];
        for (int i = 0; i < ins_cnt; ++ i) {
            points[i].obj_id = obj_id;
            points[i].ins_id = i;
            points[i].coord = new double[Dim];
            for (int j = 0; j < Dim; ++ j) {
                points[i].coord[j] = rand()/double(RAND_MAX)*(r[j] - l[j]) + l[j];
            }
            points[i].prob = 1/double(ins_cnt);
        }
        delete[] center;
        delete[] l;
        delete[] r;
        return ins_cnt;
    }

    /* generate correlated instances with unoform distribution */
    int GenCorrPoints(const int& obj_id, Instance* &points, const double& var1, const double& var2) const {
        double c = 0;
        do { c = RandNormal(0.5, var1); } while (c >= 1 || c <= 0);
        double *center = new double[Dim];
        do {
            center[Dim - 1] = Dim*c;
            for (int i = 0; i < Dim - 1; ++ i) {
                center[i] = RandNormal(0.5, var2);
                center[Dim - 1] -= center[i];
            }
        } while (center[Dim - 1] < 0 || center[Dim - 1] > 1);
        int ins_cnt = rand()%cnt + 1;
        double *l = new double[Dim];
        double *r = new double[Dim];
        for (int i = 0; i < Dim; ++ i) {
            double length = RandNormal(1.0, 0.025);
            l[i] = center[i] - length/2 > 0 ? center[i] - length/2 : 0;
            r[i] = center[i] + length/2 < 1 ? center[i] + length/2 : 1;
        }
        points = new Instance[ins_cnt];
        for (int i = 0; i < ins_cnt; ++ i) {
            points[i].obj_id = obj_id;
            points[i].ins_id = i;
            points[i].coord = new double[Dim];
            for (int j = 0; j < Dim; ++ j) {
                points[i].coord[j] = rand()/double(RAND_MAX)*(r[j] - l[j]) + l[j];
            }
            points[i].prob = 1/double(ins_cnt);
        }
        delete[] center;
        delete[] l;
        delete[] r;
        return ins_cnt;
    }
};

int main(int argc, char const *argv[]) {
    DataOperator data_operator(10, 10);
    // data_operator.GenData(1);
    data_operator.LoadData(argv[1]);
    data_operator.PrintData();
    // data_operator.WriteData(argv[1]);
    return 0;
}


#endif