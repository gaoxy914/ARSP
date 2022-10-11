#include "object.hpp"
#include "kdtree.hpp"

class Dataset {
    int dim, m, n;
    vector<int> cnt;
    vector<InstanceBase> instances;

public:
    Dataset(const int& dim, const int& m) : dim(dim), m(m), n(0) {
        cnt.resize(m);
    }

    void LoadData(const char* path) {
        ifstream file((string(path) + to_string(dim) + string("/") + to_string(m) + string("/cnt.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) { file >> cnt[i]; n += cnt[i]; }
        file.close();
        instances.resize(n);
        file.open((string(path) + to_string(dim) + string("/") + to_string(m) + string("/instances.data")).c_str(), ios::in);
        for (int i = 0; i < n; ++ i) {
            instances[i].dim = dim;
            file >> instances[i].obj_id;
            file >> instances[i].ins_id;
            file >> instances[i].prob;
            instances[i].coord = new double[dim];
            for (int j = 0; j < dim; ++ j) file >> instances[i].coord[j];
        }
        file.close();
    }

    void PrintData() {
        cout << "m = " << m << "\t n = " << n << endl;
        int s = 0;
        for (int i = 0; i < m; ++ i) {
            cout << "obj " << i << "\t cnt = " << cnt[i] << endl;
            for (int j = s; j < s + cnt[i]; ++ j) {
                cout << "id = " << instances[j].ins_id << "\t (";
                for (int k = 0; k < dim - 1; ++ k) cout << instances[j].coord[k] << ", ";
                cout << instances[j].coord[dim - 1] << ")\t";
                cout << "prob = " << 1/instances[j].prob << endl;
            }
        }
    }

    void TransAlg(const HyperBox& R, map<int, double>& results) {
        vector<vector<double>> vertices;
        int new_dim = int(pow(2, dim));
        vertices.resize(new_dim);
        for (int i = 0; i < new_dim; ++ i) {
            vertices[i].resize(dim);
            for (int j = 0; j < dim - 1; ++ j) vertices[i][j] = ((i>>j)&1) == 0 ? R.left_bottom[j] : R.right_top[j];
            vertices[i][dim - 1] = 1;
        }
        vector<InstanceBase> new_instances;
        new_instances.resize(n);
        for (int i = 0; i < n; ++ i) {
            new_instances[i].dim = new_dim;
            new_instances[i].obj_id = instances[i].obj_id;
            new_instances[i].ins_id = instances[i].ins_id;
            new_instances[i].prob = instances[i].prob;
            new_instances[i].coord = new double[new_dim];
            for (int j = 0; j < new_dim; ++ j) new_instances[i].coord[j] = instances[i].Score(vertices[j]);
        }
        // build kd-tree
        KDTree kdtree(new_dim, m, new_instances);
        // cout << "built\n";
        // calculate skyline probability
        results = kdtree.CalSkyPorb();
        // cout << "done\n";
        for (auto iter : results)
            if (iter.second != 0)
                cout << "(" << iter.first << ", " << iter.second << ")\n";
    }
};

/*
 * ./trans path dim m
 * dim: dimentionality
 * path: file path
 * m: number of uncertain tuples
 */
int main(int argc, char const *argv[]) {
    map<int, double> results;
    int dim = atoi(argv[2]);
    int m = atoi(argv[3]);
    Dataset D(dim, m);
    double *l = new double[dim - 1], *r = new double[dim - 1];
    l[0] = 1; l[1] = 1; r[0] = 2; r[1] = 2;
    HyperBox R(dim - 1, l, r);
    D.LoadData(argv[1]);
    // D.PrintData();
#ifdef _LINUX_
    struct timeval start, end;
    gettimeofday(&start, NULL);
#endif
    D.TransAlg(R, results);
#ifdef _LINUX_
    long long mtime, seconds, useconds;
    gettimeofday(&end, NULL);
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = seconds*1000000 + useconds;
    printf("Total time is: %lld\n", mtime);
#endif
    delete[] l;
    delete[] r;
    return 0;
}
