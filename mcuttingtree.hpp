#include <iostream>
#include <vector>
#include <numeric>
#include <math.h>
#include <stdlib.h>
#include <time.h>

extern "C" {
    #include "glpk.h"
}

using namespace std;

#define Dim 3

class CuttingTree {
public:
    int c;
    int n_child;
    int depth;
    double *left_bottom;
    double *right_top;

    struct Node {
        int *pivots;
        Node **children;
        vector<int> throughs;

        Node() : pivots(nullptr), children(nullptr) {}        
    };

    Node *root;
    vector<double*> planes;

    CuttingTree(int c, int n_planes, double* l, double* r) {
        srand((unsigned)time(nullptr));
        this->c = c;
        this->n_child = pow(2, c);
        this->depth = 0;
        this->left_bottom = new double[Dim];
        this->right_top = new double[Dim];
        for (int i = 0; i < Dim; ++ i) {
            this->left_bottom[i] = l[i];
            this->right_top[i] = r[i];
        }
        vector<int> through;
        through.reserve(n_planes);
        planes.reserve(n_planes);
        for (int i = 0; i < n_planes; ++ i) {
            double *p = new double[Dim];
            for (int j = 0; j < Dim; ++ j) p[j] = rand()/double(RAND_MAX);
            p[Dim - 1] *= -1;
            planes.push_back(p);
            through.push_back(i);
        }
        vector<pair<int, int> > ancestors;
        root = BuildRecursive(through, ancestors, 1);
        // cout << depth << endl;
    }

    Node* BuildRecursive(vector<int> through, vector<pair<int, int> > ancestors, int depth) {
        Node *node = new Node();
        if (depth > this->depth) this->depth = depth;
        if (through.size() <= 2*c) {
            for (auto p : through) node->throughs.push_back(p);
        } else {
            node->pivots = new int[c];
            node->children = new Node*[n_child];
            random_shuffle(through.begin(), through.end());
            for (int i = 0; i < c; ++ i) {
                node->pivots[i] = through[i];
                ancestors.push_back(make_pair(through[i], 0));
                // cout << ancestors.back().first << endl;
            }
            for (int i = 0; i < n_child; ++ i) {
                vector<int> subthrough;
                for (int j = 0; j < c; ++ j) ancestors[ancestors.size() - 1 - j].second = (i>>j)&1;
                // for (int j = ancestors.size() - c; j < ancestors.size(); ++ j) cout << ancestors[j].first << '\t' << ancestors[j].second << endl;
                /* for (int j = c; j < through.size(); ++ j) {
                    int res = Intersect(through[j], ancestors);
                    if (res == -1) { cout << "unfeasible region.\n"; break; }
                    else if (res == 1) subthrough.push_back(through[j]);
                } */
                GetIntersect(through, ancestors, subthrough);
                if (subthrough.size() > 0) {
                    // cout << "child " << i << " depth = " << depth << '\t' << "|planes| = " << subthrough.size() << endl;
                    node->children[i] = BuildRecursive(subthrough, ancestors, depth + 1);
                }
                else node->children[i] = nullptr;
            }
        }
        return node;
    }

    int Intersect(int p, vector<pair<int, int> > ancestors) {
        glp_prob *lp;
        glp_term_out(GLP_OFF);
        lp = glp_create_prob();
        glp_add_cols(lp, Dim);
        for (int i = 0; i < Dim; ++ i)
            glp_set_col_bnds(lp, i + 1, GLP_DB, left_bottom[i], right_top[i]);
        glp_set_obj_coef(lp, 1, 1);
        for (int i = 1; i < Dim; ++ i)
            glp_set_obj_coef(lp, i + 1, -planes[p][i - 1]);
        int cnt_cons = ancestors.size();
        int *ia = nullptr, *ja = nullptr;
        double *ar = nullptr;
        if (cnt_cons != 0) {
            glp_add_rows(lp, cnt_cons);
            for (int i = 0; i < cnt_cons; ++ i) {
                if (ancestors[i].second == 0) glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, planes[ancestors[i].first][Dim - 1]);
                else glp_set_row_bnds(lp, i + 1, GLP_LO, planes[ancestors[i].first][Dim - 1], 0.0);
            }
            int size = cnt_cons*Dim;
            ia = new int[size + 1];
            ja = new int[size + 1];
            ar = new double[size + 1];
            for (int i = 0; i < cnt_cons; ++ i) {
                for (int j = 0; j < Dim; ++ j) {
                    int k = i*Dim + j + 1;
                    ia[k] = i + 1;
                    ja[k] = j + 1;
                    if (j == 0) ar[k] = 1;
                    else ar[k] = -planes[ancestors[i].first][j - 1];
                }
            }
            glp_load_matrix(lp, size, ia, ja, ar);
        }
        glp_set_obj_dir(lp, GLP_MIN);
        glp_simplex(lp, NULL);
        double min_value = glp_get_obj_val(lp);
        glp_set_obj_dir(lp, GLP_MAX);
        glp_simplex(lp, NULL);
        double max_value = glp_get_obj_val(lp);
        int status = glp_get_status(lp);
        if (ia) delete[] ia;
        if (ja) delete[] ja;
        if (ar) delete[] ar;
        glp_delete_prob(lp);
        if (status == GLP_NOFEAS) return -1;
        if (max_value >= planes[p][Dim - 1] && min_value <= planes[p][Dim - 1]) return 1;
        else return 0;
    }

    void GetIntersect(vector<int> through, vector<pair<int, int> > ancestors, vector<int>& subthrough) {
        glp_prob *lp;
        glp_term_out(GLP_OFF);
        lp = glp_create_prob();
        glp_add_cols(lp, Dim);
        for (int i = 0; i < Dim; ++ i)
            glp_set_col_bnds(lp, i + 1, GLP_DB, left_bottom[i], right_top[i]);
        int cnt_cons = ancestors.size();
        glp_add_rows(lp, cnt_cons);
        for (int i = 0; i < cnt_cons; ++ i) {
            if (ancestors[i].second == 0) glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, planes[ancestors[i].first][Dim - 1]);
            else glp_set_row_bnds(lp, i + 1, GLP_LO, planes[ancestors[i].first][Dim - 1], 0.0);
        }
        int size = cnt_cons*Dim;
        int ia[size + 1];
        int ja[size + 1];
        double ar[size + 1];
        for (int i = 0; i < cnt_cons; ++ i) {
            for (int j = 0; j < Dim; ++ j) {
                int k = i*Dim + j + 1;
                ia[k] = i + 1;
                ja[k] = j + 1;
                if (j == 0) ar[k] = 1;
                else ar[k] = -planes[ancestors[i].first][j - 1];
            }
        }
        glp_load_matrix(lp, size, ia, ja, ar);
        glp_set_obj_coef(lp, 1, 1);
        for (int i = c; i < through.size(); ++ i) {
            for (int j = 1; j < Dim; ++ j)
                glp_set_obj_coef(lp, j + 1, -planes[through[i]][j - 1]);
            glp_set_obj_dir(lp, GLP_MIN);
            glp_simplex(lp, NULL);
            if (glp_get_status(lp) == GLP_NOFEAS) break;
            double min_value = glp_get_obj_val(lp);
            glp_set_obj_dir(lp, GLP_MAX);
            glp_simplex(lp, NULL);
            double max_value = glp_get_obj_val(lp);
            if (max_value >= planes[through[i]][Dim - 1] && min_value <= planes[through[i]][Dim - 1]) subthrough.push_back(through[i]);
        }
        glp_delete_prob(lp);
    }

    /* void Split(vector<int> through, vector<vector<int> >& subthroughs, vector<pair<int, int> > ancestors) {
        glp_prob *lp;
        glp_term_out(GLP_OFF);
        lp = glp_create_prob();
        glp_add_cols(lp, Dim);
        for (int i = 0; i < Dim; ++ i)
            glp_set_col_bnds(lp, i + 1, GLP_DB, left_bottom[i], right_top[i]);
        glp_set_obj_coef(lp, 1, 1);
    } */
};

int main(int argc, char const *argv[]) {
    double *l = new double[Dim], *r = new double[Dim];
    for (int i = 0; i < Dim; ++ i) { l[i] = -10; r[i] = 0; }
    int n_planes = 1000000;
    CuttingTree tree(Dim*Dim, n_planes, l, r);
    cout << tree.depth << endl;
    /* double c1[2] = {-1, 1.5};
    double c2[2] = {-0.5, 0.7};
    vector<pair<double*, int> > ancestors;
    ancestors.push_back(make_pair(c1, 0));
    ancestors.push_back(make_pair(c2, 1));
    double p[2] = {-1, 1.6};
    cout << Intersect(l, r, p, ancestors) << endl; */
    return 0;
}