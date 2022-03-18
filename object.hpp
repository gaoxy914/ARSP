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

#define Dim 3

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
    int d;
    double *left_bottom, *right_top;

    HyperBox() {
        left_bottom = right_top = nullptr;
    }

    HyperBox(const int& d, double* coord) {
        this->d = d;
        this->left_bottom = new double[d];
        this->right_top = new double[d];
        for (int i = 0; i < d; ++ i) {
            this->left_bottom[i] = coord[i];
            this->right_top[i] = coord[i];
        }
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

    /* return centroid of the box */
    double* GetCentroid() const {
        double *centroid = new double[d];
        for (int i = 0; i < d; ++ i) {
            centroid[i] = (left_bottom[i] + right_top[i])/2;
        }
        return centroid;
    }

    /* return 2^d extremes */
    void GetExtremes(vector<double*>& extremes) const {
        extremes.resize(int(pow(2, d)));
        for (int i = 0; i < pow(2, d); ++ i) {
            int k = i;
            extremes[i] = new double[d + 1];
            for (int j = 0; j < d; ++ j) {
                extremes[i][j] = ((k>>j)&1) == 0 ? left_bottom[j] : right_top[j];
            }
            extremes[i][d] = 1;
        }
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

    /* true : right_top R-dominates instance */
    bool RDominates(const double* coord, const HyperBox& box) const {
        double sum = coord[Dim - 1] - right_top[Dim - 1];
        for (int i = 0; i < Dim - 1; ++ i) {
            if (coord[i] - right_top[i] > 0) sum += (coord[i] - right_top[i])*box.left_bottom[i];
            else sum += (coord[i] - right_top[i])*box.right_top[i];
        }
        return sum >= 0;
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

    /* dominate test for loop based algorithm */
    bool Dominates(const InstanceBase& instance) const {
        if (ins_id == instance.ins_id) return false;
        for (int i = 0; i < Dim; ++ i) if (coord[i] > instance.coord[i]) return false;
        return true;
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

    /* R-dominate test */
    bool RDominates(const InstanceBase& instance, const HyperBox& box) const {
        double sum = instance.coord[Dim - 1] - coord[Dim - 1];
        for (int i = 0; i < Dim - 1; ++ i) {
            if (instance.coord[i] - coord[i] > 0) sum += (instance.coord[i] - coord[i])*box.left_bottom[i];
            else sum += (instance.coord[i] - coord[i])*box.right_top[i];
        }
        return sum >= 0;
    }
};


class DataOperator {
    int m, n, cnt;
    int* ins_cnt;
    vector<InstanceBase*> dataset;
public:
    DataOperator() : m(0), n(0), cnt(0) {}

    DataOperator(const int& m) : m(m), n(0), cnt(0) {}

    DataOperator(const int& m, const int& cnt) {
        this->m = m;
        this->n = 0;
        this->cnt = cnt;
        this->ins_cnt = new int[m];
        dataset.resize(m);
    }

    void PrintData() {
        cout << "m = " << m << "\t n = " << n << endl;
        for (int i = 0; i < m; ++ i) {
            cout << "obj " << i << "\t cnt = " << ins_cnt[i] << endl;
            for (int j = 0; j < ins_cnt[i]; ++ j) {
                cout << "id = " << dataset[i][j].ins_id << "\t (";
                for (int k = 0; k < Dim - 1; ++ k) cout << dataset[i][j].coord[k] << ", ";
                cout << dataset[i][j].coord[Dim - 1] << ")\t";
                cout << "prob = " << 1/dataset[i][j].prob << endl;
            }
        }
    }

    void LoadData(const char *path) {
        ifstream file((string(path) + to_string(m) + string("/cnt.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) { file >> ins_cnt[i]; n += ins_cnt[i]; }
        file.close();
        file.open((string(path) + to_string(m) + string("/instances.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) {
            dataset[i] = new InstanceBase[ins_cnt[i]];
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
            for (int i = 0; i < m; ++ i)
                ins_cnt[i] = GenIndePoints(i, dataset[i]);
            break;
        case 2:
            for (int i = 0; i < m; ++ i)
                ins_cnt[i] = GenAntiPoints(i, dataset[i], var1);
            break;
        case 3:
            for (int i = 0; i < m; ++ i)
                ins_cnt[i] = GenCorrPoints(i, dataset[i], var1, var2);
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
    int GenIndePoints(const int& obj_id, InstanceBase* &points) {
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
        points = new InstanceBase[ins_cnt];
        for (int i = 0; i < ins_cnt; ++ i) {
            points[i].obj_id = obj_id;
            // points[i].ins_id = i;
            points[i].ins_id = n ++;
            points[i].coord = new double[Dim];
            for (int j = 0; j < Dim; ++ j) points[i].coord[j] = rand()/double (RAND_MAX)*(r[j] - l[j]) + l[j];
            // points[i].prob = 1/double(ins_cnt);
            points[i].prob = double(ins_cnt);
        }
        delete[] center;
        delete[] l;
        delete[] r;
        return ins_cnt;
    }

    /* generate anti-correlated instances with uniform distribution */
    int GenAntiPoints(const int& obj_id, InstanceBase* &points, const double& var) {
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
        points = new InstanceBase[ins_cnt];
        for (int i = 0; i < ins_cnt; ++ i) {
            points[i].obj_id = obj_id;
            // points[i].ins_id = i;
            points[i].ins_id = n ++;
            points[i].coord = new double[Dim];
            for (int j = 0; j < Dim; ++ j) {
                points[i].coord[j] = rand()/double(RAND_MAX)*(r[j] - l[j]) + l[j];
            }
            // points[i].prob = 1/double(ins_cnt);
            points[i].prob = double(ins_cnt);
        }
        delete[] center;
        delete[] l;
        delete[] r;
        return ins_cnt;
    }

    /* generate correlated instances with unoform distribution */
    int GenCorrPoints(const int& obj_id, InstanceBase* &points, const double& var1, const double& var2) {
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
        points = new InstanceBase[ins_cnt];
        for (int i = 0; i < ins_cnt; ++ i) {
            points[i].obj_id = obj_id;
            // points[i].ins_id = i;
            points[i].ins_id = n ++;
            points[i].coord = new double[Dim];
            for (int j = 0; j < Dim; ++ j) {
                points[i].coord[j] = rand()/double(RAND_MAX)*(r[j] - l[j]) + l[j];
            }
            // points[i].prob = 1/double(ins_cnt);
            points[i].prob = double(ins_cnt);
        }
        delete[] center;
        delete[] l;
        delete[] r;
        return ins_cnt;
    }
};

/* int main(int argc, char const *argv[]) {
    DataOperator data_operator(10, 10);
    data_operator.GenData(1);
    // data_operator.LoadData(argv[1]);
    data_operator.PrintData();
    data_operator.WriteData(argv[1]);
    return 0;
} */


#endif