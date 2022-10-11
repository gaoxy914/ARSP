#include "object.hpp"

class Instance : public InstanceBase {
public:
    double beta;

    Instance() : InstanceBase(), beta(1) {}

    Instance& operator =(const Instance& instance) {
        if (&instance != this) {
            InstanceBase::operator=(instance);
            beta = instance.beta;
        }
        return *this;
    }

    double GetProb() const { return (1/prob)*beta; }
};

struct InstanceComparator {
    InstanceComparator(double *_weight) : weight(_weight) {}

    bool operator ()(const Instance& u, const Instance& v) const {
        return u.Score(weight) < v.Score(weight);
    }
    double *weight;
};

class Dataset {
    int dim, m, n, c;
    int* cnt;
    vector<Instance> instances;
    map<pair<int, int>, bool> dominate;

public:
    Dataset(const int& dim, const int& m, const int& c) : n(0), dim(dim), m(m), c(c) {
        this->cnt = new int[m];
        memset(this->cnt, 0, m*sizeof(int));
    }

    void LoadData(const char* path) {
        ifstream file((string(path) + to_string(dim) + string("/") + to_string(m) + string("/cnt.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) { file >> cnt[i]; n += cnt[i]; }
        file.close();
        file.open((string(path) + to_string(dim) + string("/") + to_string(m) + string("/instances.data")).c_str(), ios::in);
        instances.resize(n);
        for (int i = 0; i < n; ++ i) {
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
            s += cnt[i];
        }
    }

    void LoopPreprocessing() {
        double* weight = new double[dim];
        memset(weight, 0, dim*sizeof(double));
        weight[0] = 1;
        sort(instances.begin(), instances.end(), InstanceComparator(weight));
        int* sigma = new int[m];
        for (int i = 0; i < instances.size(); ++ i) {
            memset(sigma, 0, m*sizeof(int));
            vector<int> dominators;
            for (int j = 0; j < i - 1; ++ j) {
                if (instances[i].obj_id == instances[j].obj_id) continue;
                if (instances[j].Dominates(instances[i])) {
                    if (instances[i].ins_id == 60) cout << instances[j].ins_id << endl;
                    dominators.push_back(instances[j].ins_id);
                    sigma[instances[j].obj_id] ++;
                    if (sigma[instances[j].obj_id] == cnt[instances[j].obj_id]) {
                        instances[i].beta = 0;
                        break;
                    }
                }
            }
            if (instances[i].beta != 0) {
                for (int id : dominators) {
                    pair<int, int> p = make_pair(id, instances[i].ins_id);
                    dominate[p] = true;
                }
            }
        }
        delete[] sigma;
        delete[] weight;
    }

    void LoopAlg(const HyperBox& R, map<int, double>& results) {
        // vector<double*> MBB; // minimum bounding rectangle with capacity c
        double* weight = new double[dim];
        for (int i = 0; i < dim - 1; ++ i) weight[i] = (R.left_bottom[i] + R.right_top[i])/2;
        weight[dim - 1] = 1;
        sort(instances.begin(), instances.end(), InstanceComparator(weight));
        int* sigma = new int[m];
        for (int i = 0; i < n; ++ i) {
            memset(sigma, 0, m*sizeof(int));
            double beta = instances[i].beta;
            for (int j = 0; j < n; ++ j) {
                if (instances[j].obj_id != instances[i].obj_id) {
                    if (instances[j].RDominates(instances[i], R)) {
                        beta *= (instances[j].prob - sigma[instances[j].obj_id] - 1)/(instances[j].prob - sigma[instances[j].obj_id]);
                        sigma[instances[j].obj_id] ++;
                    }
                }
            }
            cout << endl;
            results[instances[i].ins_id] = (1/instances[i].prob)*beta;
        }
        delete[] weight;
        delete[] sigma;
        for (auto iter : results)
            if (iter.second != 0)
                cout << "(" << iter.first << ", " << iter.second << ")\n";
    }
};

/*
 * ./loop path dim m c
 * path: file path
 * m: number of uncertain tuples
 * c: capacity of MBB
 */
int main(int argc, char const *argv[]) {
    map<int, double> results;
    int dim = atoi(argv[2]);
    int m = atoi(argv[3]);
    int c = atoi(argv[4]);
    Dataset D(dim, m, c);
    // double l[dim - 1] = {1, 1};
    // double r[dim - 1] = {2, 2};
    // HyperBox R(dim - 1, l, r);
    D.LoadData(argv[1]);
    D.PrintData();
#ifdef _LINUX_
    struct timeval start, end1, end2;
    gettimeofday(&start, NULL);
#endif
    D.LoopPreprocessing();
#ifdef _LINUX_
    gettimeofday(&end1, NULL);
    long long mtime, seconds, useconds;
    seconds = end1.tv_sec - start.tv_sec;
    useconds = end1.tv_usec - start.tv_usec;
    mtime = seconds*1000000 + useconds;
    printf("Preprocessing Time is : %lld\n", mtime);
#endif
    // D.LoopAlg(R, results);
#ifdef _LINUX_
    gettimeofday(&end2, NULL);
    seconds = end2.tv_sec - start.tv_sec;
    useconds = end2.tv_usec - start.tv_usec;
    mtime = seconds*1000000 + useconds;
    printf("Total time is: %lld\n", mtime);
#endif
    return 0;
}


