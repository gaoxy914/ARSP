#include <iostream>
#include <vector>
#include <numeric>

extern "C" {
    #include "glpk.h"
}

using namespace std;

#define Dim 2

//对于d的扩展性很差，d越大，越难通过本身的平面来将它们划分开

// x[d] = w[1]x[1] + ... + w[d-1]x[d-1] + w[d]
class HyperPlane {
public:
    int id;
    double *coef; // w[1] ... w[d]

    HyperPlane() : id(-1), coef(nullptr) {}

};

class CuttingTree {
public:
    int c;
    int depth;
    double *left_bottom;
    double *right_top;

    struct Node {
        int id;
        vector<int> throughs;
        Node *left, *right; // left: planes below id, right: planes above id

        Node() : id(-1), left(nullptr), right(nullptr) {}
    };

    Node *root;    
    vector<HyperPlane> planes;

    CuttingTree(const int& c, const int& n_planes, const double* l, const double* r) {
        this->c = c;
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
            if (i%10000 == 0) cout << i << endl;
            HyperPlane p;
            p.id = i;
            p.coef = new double[Dim];
            for (int j = 0; j < Dim; ++ j) p.coef[j] = rand()/double(RAND_MAX);
            p.coef[Dim - 1] *= -1;
            planes.push_back(p);
            through.push_back(i);
        }
        vector<pair<int, int> > ancestors;
        root = BuildRecursive(through, ancestors, 1);
    }

    Node* BuildRecursive(vector<int> through, vector<pair<int, int> > ancestors, int depth) {
        // cout << "|through| = " << through.size() << endl;
        Node *node = new Node();
        if (depth > this->depth) this->depth = depth;
        if (through.size() <= c) {
            for (auto p : through) node->throughs.push_back(p);
        } else {
            vector<int> left, right;
            node->id = __choose(through, ancestors, left, right);
            // cout << "finding split plane.\n";
            // node->id = __rand_choose(through, ancestors, left, right);
            // cout << "split at plane " << node->id << endl;
            cout << "depth = " << depth << '\t' << "left = " << left.size() << '\t' << "right = " << right.size() << endl;
            ancestors.push_back(make_pair(node->id, -1));
            node->left = BuildRecursive(left, ancestors, depth + 1);
            ancestors.back().second = 1;
            node->right = BuildRecursive(right, ancestors, depth + 1);
            ancestors.pop_back();
        }
        return node;
    }

    int __choose(const vector<int>& through, const vector<pair<int, int> > ancestors, vector<int>& left, vector<int>& right) const {
        int p = -1, min_sum = through.size();
        for (int i = 0; i < through.size(); ++ i) {
            int sum = 0;
            vector<int> l, r;
            for (int j = 0; j < through.size(); ++ j) {
                if (i == j) continue;
                ++ sum;
                int relation = __relation(through[i], through[j], ancestors);
                if (relation == 0) ++ sum;
                if (relation < 1) r.push_back(through[j]);
                if (relation > -1) l.push_back(through[j]);
            }
            if (max(l.size(), r.size()) < min_sum) {
                p = through[i];
                min_sum = max(l.size(), r.size());
                left.swap(l);
                right.swap(r);
            }
        }
        return p;
    }

    int __rand_choose(const vector<int>& through, const vector<pair<int, int> > ancestors, vector<int>& left, vector<int>& right) const {
        int p = -1, min_sum = through.size();
        for (int i = 0; i < 20; ++ i) {
            // cout << "i = " << i << endl;
            int s = rand()%through.size();
            // cout << "random s = " << s << endl;
            int sum = 0;
            vector<int> l, r;
            for (int j = 0; j < through.size(); ++ j) {
                if (s == j) continue;
                ++ sum;
                // cout << s << '\t' << j << endl;
                // cout << "compare " << through[s] << " and " << through[j] << endl;
                int relation = __relation(through[s], through[j], ancestors);
                // cout << "done.\n";
                // cout << relation << endl;
                if (relation == 0) ++ sum;
                if (relation < 1) r.push_back(through[j]);
                if (relation > -1) l.push_back(through[j]);
            }
            if (max(l.size(), r.size()) < min_sum) {
                p = through[s];
                min_sum = max(l.size(), r.size());
                left.swap(l);
                right.swap(r);
            }
        }
        return p;
    }

    /*
     * constrcut LP problem with constraints induced by ancestors, and the object function if p - q
     * return : 1 if min >= 0 (p above q), -1 else if max <= 0(p below q), 0 otherwise (p intersects q)
     */
    int __relation(const int& p, const int& q, vector<pair<int, int> > ancestors) const {
        glp_prob *lp;
        glp_term_out(GLP_OFF);
        lp = glp_create_prob();
        glp_add_cols(lp, Dim);
        for (int i = 0; i < Dim; ++ i)
            glp_set_col_bnds(lp, i + 1, GLP_DB, left_bottom[i], right_top[i]);
        glp_set_obj_coef(lp, 1, 0);
        for (int i = 1; i < Dim; ++ i)
            glp_set_obj_coef(lp, i + 1, planes[p].coef[i - 1] - planes[q].coef[i - 1]);
        int cnt_cons = ancestors.size();
        int *ia = nullptr, *ja = nullptr;
        double *ar = nullptr;
        if (cnt_cons != 0) {
            glp_add_rows(lp, cnt_cons);
            for (int i = 0; i < cnt_cons; ++ i) {
                if (ancestors[i].second == 0) glp_set_row_bnds(lp, i + 1, GLP_UP, 0.0, planes[ancestors[i].first].coef[Dim - 1]);
                else glp_set_row_bnds(lp, i + 1, GLP_LO, planes[ancestors[i].first].coef[Dim - 1], 0.0);
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
                    else ar[k] = -planes[ancestors[i].first].coef[j - 1];
                }
            }
            glp_load_matrix(lp, size, ia, ja, ar);
        }
        glp_set_obj_dir(lp, GLP_MIN);
        glp_simplex(lp, NULL);
        double min_value = glp_get_obj_val(lp);
        min_value += (planes[p].coef[Dim - 1] - planes[q].coef[Dim - 1]);
        glp_set_obj_dir(lp, GLP_MAX);
        glp_simplex(lp, NULL);
        double max_value = glp_get_obj_val(lp);
        max_value += (planes[p].coef[Dim - 1] - planes[q].coef[Dim - 1]);
        if (ia) delete[] ia;
        if (ja) delete[] ja;
        if (ar) delete[] ar;
        glp_delete_prob(lp);
        if (min_value >= 0) return 1;
        else if (max_value <= 0) return -1;
        else return 0;
    }
};

int main(int argc, char const *argv[]) {
    double *l = new double[Dim], *r = new double[Dim];
    for (int i = 0; i < Dim; ++ i) { l[i] = -10; r[i] = 0; }
    int n_planes = 1000;
    CuttingTree tree(50, n_planes, l, r);
    cout << tree.depth << endl;
    return 0;
}
