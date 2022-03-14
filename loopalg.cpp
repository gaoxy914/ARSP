#include "object.cpp"

struct InstanceComparator {
    InstanceComparator(double *_weight) : weight(_weight) {}
    bool operator ()(const Instance& u, const Instance& v) const {
        return dot(u.coord, weight) < dot(v.coord, weight);
    }
    double *weight;
};

class Instance : public InstanceBase {
public:
    double beta;
    double *sigma;

    Instance() : InstanceBase(), beta(1), sigma(nullptr) {}

    virtual ~Instance() {
        if (sigma != nullptr) delete[] sigma;
        sigma = nullptr;
    }

    void InitialSigma(const int& m) {
        sigma = new double[m];
        for (int i = 0; i < m; ++ i) sigma[i] = 1;
    }
};

class Object {
public:
    int obj_id;
    int cnt;
    double eprob;
    Instance* instances;

    vector<Instance> heap; // for local sort

    Object(const int& obj_id) : obj_id(obj_id), eprob(0), instances(nullptr) {}

    ~Object() {
        if (instances != nullptr) delete[] instances;
        instances = nullptr;
    }

    void sort(double *weight) {
        heap.reserve(cnt);
        for (int i = 0; i < cnt; ++ i) heap.push_back(instances[i]);
        make_heap(heap.begin(), heap.end(), InstanceComparator(weight));
    }

    Instance pop() {
        Instance instance = heap.front();
        pop_heap(heap.begin(), heap.end());
        heap.pop_back();
        return instance;
    }

    Instance top() {
        return heap.front();
    }
};

class Dataset {
    int m, n;
    vector<Object> objects;
    map<pair<int, int>, bool> dominate;

public:
    Dataset(const int& m) {
        n = 0;
        this->m = m;
        objects.resize(n);
    }

    void LoadData(const char* path) {
        ifstream file((string(path) + to_string(m) + string("/cnt.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) { file >> objects[i].cnt; n += objects[i].cnt; }
        file.close();
        file.open((string(path) + to_string(m) + string("/instances.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) {
            objects[i].instances = new Instance[objects[i].cnt];
            for (int j = 0; j < objects[i].cnt; ++ j) {
                file >> objects[i].instances[j].obj_id;
                file >> objects[i].instances[j].ins_id;
                file >> objects[i].instances[j].prob;
                objects[i].instances[j].coord = new double[Dim];
                for (int k = 0; k < Dim; ++ k) file >> objects[i].instances[j].coord[k];
            }
        }
        file.close();        
    }

    void LoopPreprocessing() {
        for (int i = 0; i < m - 1; ++ i) {
            Object obj1 = objects[i];
            for (int p = 0; p < obj1.cnt; ++ p) {
                for (int j = i + 1; j < m; ++ j) {
                    Object obj2 = objects[j];
                    for (int q = 0; q < obj2.cnt; ++ q) {
                        if (obj1.instances[p].Dominates(obj2.instances[q]))
                            dominate[make_pair(obj1.instances[p].ins_id, obj2.instances[q].ins_id)] = true;
                        else if (obj2.instances[q].Dominates(obj1.instances[p]))
                            dominate[make_pair(obj2.instances[q].ins_id, obj1.instances[p].ins_id)] = true;
                    }
                }
            }
        }
    }

    void LoopAlg(const HyperBox& R) {
        vector<Instance> min_heap; // heap for global sorting
        vector<Instance> popped_instances; // instances popped before the current one
        min_heap.reserve(m);
        double* centroid = R.GetCentroid();
        double* weight = new double[Dim];
        for (int i = 0; i < Dim - 1; ++ i) weight[i] = centroid[i];
        weight[Dim - 1] = 1;
        for (int i = 0; i < n; ++ i) {
            objects[i].sort(weight);
            Instance ins = objects[i].pop();
            min_heap.push_back(ins);
        }
        make_heap(min_heap.begin(), min_heap.end(), InstanceComparator(weight));        
        while (min_heap.size() > 0) {
            Instance cur_ins = min_heap.front();
            pop_heap(min_heap.begin(), min_heap.end());
            min_heap.pop_back();

        }
        delete[] centroid;
        delete[] weight;
    }
};

int main(int argc, char const *argv[]) {
    // Dataset dataset();
    // dataset.gen
#ifdef _LINUX_
    struct timeval start, end;
    gettimeofday(&start, NULL);
#endif
    // dataset.loop_bl();
#ifdef _LINUX_
    gettimeofday(&end, NULL);
    long long mtime, seconds, useconds;
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = seconds*1000000 + useconds;
    printf("Total time is: %lld\n", mtime);
#endif
    return 0;
}


