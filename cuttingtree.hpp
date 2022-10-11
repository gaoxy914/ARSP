#include "object.hpp"

class CuttingTree {
    int c;
    int n_child;
    int depth;
    HyperBox R;

    struct Node {
        int *pivots;
        Node **children;
        int n_planes;
        int *planes;

        Node() : pivots(nullptr), children(nullptr), n_planes(0), planes(nullptr) {}
    };

    Node *root;
    vector<HyperPlane> planes;

    struct Info {
        int* sigma;
        double beta;
        int xi;
    };

    Node* BuildRecursive(vector<int> through, vector<pair<int, int> > constraints, const int& depth) {
        Node *node = new Node();
        if (through.size() <= c) {
            node->n_planes = through.size();
            node->planes = new int[through.size()];
            for (int i = 0; i < node->n_planes; ++ i) node->planes[i] = through[i];
        } else {
            random_shuffle(through.begin(), through.end());
            node->pivots = new int[c];
            for (int i = 0; i < c; ++ i) {
                node->pivots[i] = through[i];
                constraints.push_back(make_pair(through[i], 0));
            }
            node->children = new Node*[n_child];
            for (int i = 0; i < n_child; ++ i) {
                vector<int> subthrough;
                for (int j = 0; j < c; ++ j) constraints[constraints.size() - 1 - j].second = (i>>j)&1;
                for (int j = c; j < through.size(); ++ j) if (Intersect(through[j], constraints)) subthrough.push_back(through[j]);
                node->children[i] = BuildRecursive(subthrough, constraints, depth + 1);
            }
        }
    }

    bool Intersect(const int& p, const vector<pair<int, int> >& constraints) const {
        glp_prob *lp;
        glp_term_out(GLP_OFF);
        lp = glp_create_prob();
        glp_add_cols(lp, Dim);
        for (int i = 0; i < Dim; ++ i) glp_set_col_bnds(lp, i + 1, GLP_DB, R.left_bottom[i], R.right_top[i]);
        glp_set_obj_coef(lp, 1, 1);
        for (int i = 1; i < Dim; ++ i) glp_set_obj_coef(lp, i + 1, -planes[p].coef[i - 1]);
        int n_constraints = constraints.size();
        int *ia = nullptr, *ja = nullptr;
        double *ar = nullptr;
        if (n_constraints != 0) {
            glp_add_rows(lp, n_constraints);
            for (int i = 0; i < n_constraints; ++ i) {
                if (constraints[i].second == 0) glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, planes[constraints[i].first].coef[Dim - 1]);
                else glp_set_row_bnds(lp, i + 1, GLP_LO, planes[constraints[i].first].coef[Dim - 1], 0.0);
            }
            int size = n_constraints*Dim;
            ia = new int[size + 1];
            ja = new int[size + 1];
            ar = new double[size + 1];
            for (int i = 0; i < n_constraints; ++ i) {
                for (int j = 0; j < Dim; ++ j) {
                    int k = i*Dim + j + 1;
                    ia[k] = i + 1;
                    ja[k] = j + 1;
                    if (j == 0) ar[k] = 1;
                    else ar[k] = -planes[constraints[i].first].coef[j - 1];
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
        if (ia) delete[] ia;
        if (ja) delete[] ja;
        if (ar) delete[] ar;
        glp_delete_prob(lp);
        if (max_value > planes[p].coef[Dim - 1] && min_value < planes[p].coef[Dim - 1]) return true;
        else return false;
    }

public:
    CuttingTree() : root(nullptr), c(0), n_child(0), depth(0) {}
    CuttingTree(const int& c, const HyperBox& R, const vector<HyperPlane>& planes) {
        this->root = nullptr;
        this->c = c;
        this->n_child = pow(2, c);
        this->planes = planes;
        this->R = R;
        
    }
};