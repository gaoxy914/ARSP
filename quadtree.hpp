#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>

using namespace std;

#define Dim 3

class QuadTree {
public:
    int c;
    int n_child = pow(2, Dim);
    double *left_bottom;
    double *right_top;
    int depth;
    vector<double*> planes;

    struct Node {
        Node **children;
        int n_planes;
        int *planes;
        double *centor; // try random centor

        Node() : n_planes(0), children(nullptr), planes(nullptr) {}
    };

    Node *root;

    void GetSubSpace(const double* l, const double* r, double* centor, const int& k, double*& subl, double*& subr) const {
        for (int i = 0; i < Dim; ++ i) {
            if ((k>>i)&1) {
                subl[i] = centor[i];
                subr[i] = r[i];
            } else {
                subl[i] = l[i];
                subr[i] = centor[i];
            }
            // subl[i] = l[i] + (r[i] - l[i])/2*((k>>i)&1);
            // subr[i] = subl[i] + (r[i] - l[i])/2;
        }
    }

    void GetSubSpace(const double* l, const double* r, const int& k, double*& subl, double*& subr) const {
        for (int i = 0; i < Dim; ++ i) {
            subl[i] = l[i] + (r[i] - l[i])/2*((k>>i)&1);
            subr[i] = subl[i] + (r[i] - l[i])/2;
        }
    }

    bool Intersect(const double* l, const double* r, const double* p) const {
        double max_value = r[Dim - 1];
        double min_value = l[Dim - 1];
        for (int i = 0; i < Dim - 1; ++ i) {
            max_value -= p[i]*l[i];
            min_value -= p[i]*r[i];
        }
        return max_value > p[Dim - 1] && min_value < p[Dim - 1];
    }

    Node* BuildRecursive(const double* l, const double* r, vector<int> through, int depth) {
        Node *node = new Node();
        this->depth = max(this->depth, depth);
        cout << depth << endl;
        if (through.size() <= c) {
            node->n_planes = through.size();
            if (node->planes != 0) for (int i = 0; i < through.size(); ++ i) node->planes[i] = through[i];
        } else {
            node->children = new Node*[n_child];
            node->centor = new double[Dim];
            for (int i = 0; i < Dim; ++ i) node->centor[i] = drand48()*(r[i] - l[i]) + l[i];
            for (int i = 0; i < n_child; ++ i) {
                double *subl = new double[Dim], *subr = new double[Dim];
                GetSubSpace(l, r, i, subl, subr);
                vector<int> subthrough;
                for (int p : through)
                    if (Intersect(subl, subr, planes[p])) subthrough.push_back(p);
                cout << "depth = " << depth << ", " << i << "-th children |subthrough| = " << subthrough.size() << endl;
                node->children[i] = BuildRecursive(subl, subr, subthrough, depth + 1);
                delete[] subl;
                delete[] subr;
            }
        }
        return node;
    }

    QuadTree(const int& c, const int& n_planes, const double* l, const double* r) {
        this->c = c;
        this->left_bottom = new double[Dim];
        this->right_top = new double[Dim];
        for (int i = 0; i < Dim; ++ i) {
            this->left_bottom[i] = l[i];
            this->right_top[i] = r[i];
        }
        this->depth = 0;
        vector<int> through;
        through.reserve(n_planes);
        planes.reserve(n_planes);
        for (int i = 0; i < n_planes; ++ i) {
            double* p = new double[Dim];
            for (int j = 0; j < Dim; ++ j) p[j] = rand()/double(RAND_MAX);
            p[Dim - 1] *= -1;
            planes.push_back(p);
            through.push_back(i);
        }
        root = BuildRecursive(l, r, through, 1);
        cout << depth << endl;
    }
};

int main(int argc, char const *argv[]) {
    srand((unsigned)time(nullptr));
    double *l = new double[Dim], *r = new double[Dim];
    for (int i = 0; i < Dim; ++ i) { l[i] = -10; r[i] = 0; }
    int n_planes = 1000;
    QuadTree tree(3, n_planes, l, r);
    return 0;
}
