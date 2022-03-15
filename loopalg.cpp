#include "object.cpp"

class Instance : public InstanceBase {
public:
    double beta;
    double *sigma;

    Instance() : InstanceBase(), beta(1), sigma(nullptr) {}

    void Initial(const int& m) {
        if (sigma == nullptr) sigma = new double[m];
        memset(sigma, 0, m*sizeof(sigma));
        beta = 1;
    }

    double GetProb() const {
        return prob*beta;
    }
};

struct InstanceComparator {
    InstanceComparator(double *_weight) : weight(_weight) {}
    bool operator ()(const Instance& u, const Instance& v) const {
        return dot(u.coord, weight) > dot(v.coord, weight);
    }
    double *weight;
};

class Object {
public:
    int obj_id;
    int cnt;
    int heap_size;
    double eprob;
    vector<Instance> instances;

    Object() : obj_id(-1), cnt(0), heap_size(0), eprob(0) {}

    Object(const int& obj_id) : obj_id(obj_id), cnt(0), heap_size(0), eprob(0) {}

    void MakeHeap(double *weight) {
        heap_size = cnt;
        make_heap(instances.begin(), instances.begin() + heap_size, InstanceComparator(weight));
    }

    Instance Pop(double* weight) {
        Instance instance = instances.front();
        pop_heap(instances.begin(), instances.begin() + heap_size, InstanceComparator(weight));
        heap_size --;
        return instance;
    }

    Instance Top() {
        return instances.front();
    }

    bool Empty() {
        return heap_size == 0;
    }
};

class Dataset {
    int m, n, c;
    vector<Object> objects;
    map<pair<int, int>, bool> dominate;

public:
    Dataset(const int& m, const int& c) {
        n = 0;
        this->m = m;
        this->c = c;
        objects.resize(m);
    }

    void LoadData(const char* path) {
        ifstream file((string(path) + to_string(m) + string("/cnt.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) { file >> objects[i].cnt; n += objects[i].cnt; }
        file.close();
        file.open((string(path) + to_string(m) + string("/instances.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) {
            // objects[i].instances = new Instance[objects[i].cnt];
            objects[i].instances.resize(objects[i].cnt);
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

    void PrintData() {
        cout << "m = " << m << "\t n = " << n << endl;
        for (int i = 0; i < m; ++ i) {
            cout << "obj " << i << "\t cnt = " << objects[i].cnt << endl;
            for (int j = 0; j < objects[i].cnt; ++ j) {
                cout << "id = " << objects[i].instances[j].ins_id << "\t (";
                for (int k = 0; k < Dim - 1; ++ k) cout << objects[i].instances[j].coord[k] << ", ";
                cout << objects[i].instances[j].coord[Dim - 1] << ")\t";
                cout << "prob = " << objects[i].instances[j].prob << endl;
            }
        }
    }

    void LoopPreprocessing() {
        for (int i = 0; i < m - 1; ++ i) {
            Object obj1 = objects[i];
            for (int j = i + 1; j < m; ++ j) {
                Object obj2 = objects[j];
                for (int p = 0; p < obj1.cnt; ++ p) {
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

    void LoopAlg(const HyperBox& R, map<int, double> results) {
        vector<Instance> min_heap; // heap for global sorting
        vector<Instance> popped; // instances popped before the current one
        vector<HyperBox> MBR; // minimum bounding rectangle with capacity c
        int cur_n = 0;
        min_heap.reserve(m);
        double* centroid = R.GetCentroid();
        double* weight = new double[Dim];
        for (int i = 0; i < Dim - 1; ++ i) weight[i] = centroid[i];
        weight[Dim - 1] = 1;
        for (int i = 0; i < m; ++ i) {
            objects[i].MakeHeap(weight);
            Instance ins = objects[i].Pop(weight);
            min_heap.push_back(ins);
        }
        make_heap(min_heap.begin(), min_heap.end(), InstanceComparator(weight));        
        // cout << "global heap done.\n";
        while (min_heap.size() > 0) {
            Instance cur_ins = min_heap.front();
            pop_heap(min_heap.begin(), min_heap.end(), InstanceComparator(weight));
            min_heap.pop_back();
            cur_ins.Initial(m);
            for (int i = 0; i < MBR.size(); ++ i) {
                if (MBR[i].RDominates(cur_ins.coord, R)) {
                    for (int j = i*c; j < min((i + 1)*c, cur_n); ++ j) {
                        if (popped[j].obj_id == cur_ins.obj_id) continue;
                        Instance cmp_ins = popped[j];
                        cur_ins.beta *= (1 - cur_ins.sigma[cmp_ins.obj_id] - popped[j].prob)/(1 - cur_ins.sigma[popped[j].obj_id]);
                        cur_ins.sigma[popped[j].obj_id] += popped[j].prob;
                    }
                } else {
                    for (int j = i*c; j < min((i + 1)*c, cur_n); ++ j) {
                        if (popped[j].obj_id == cur_ins.obj_id) continue;
                        if (dominate[make_pair(popped[j].ins_id, cur_ins.ins_id)] || popped[j].RDominates(cur_ins, R)) {                            
                            cur_ins.beta *= (1 - cur_ins.sigma[popped[j].obj_id] - popped[j].prob)/(1 - cur_ins.sigma[popped[j].obj_id]);
                            cur_ins.sigma[popped[j].obj_id] += popped[j].prob;
                        }
                    }
                }
            }
            /* for (auto iter : popped) { // trivial loop
                if (iter.RDominates(cur_ins, R)) {
                    cur_ins.beta *= (1 - cur_ins.sigma[iter.obj_id] - iter.prob)/(1 - cur_ins.sigma[iter.obj_id]);
                    cur_ins.sigma[iter.obj_id] += iter.prob; 
                }
            } */
            // cout << "(" << cur_ins.ins_id << ", " << cur_ins.GetProb() << ")\n";
            results[cur_ins.ins_id] = cur_ins.GetProb();
            cur_n ++;
            HyperBox B = HyperBox(Dim, cur_ins.coord, cur_ins.coord);
            if (cur_n%c == 1) MBR.push_back(B);
            else MBR.back() = MBR.back() + B;
            popped.push_back(cur_ins);
            if (!objects[cur_ins.obj_id].Empty()) {
                Instance next_ins = objects[cur_ins.obj_id].Pop(weight);
                min_heap.push_back(next_ins);
                push_heap(min_heap.begin(), min_heap.end(), InstanceComparator(weight));
            }
        }
        delete[] centroid;
        delete[] weight;
        for (auto iter : objects) {
            for (int i = 0; i < iter.cnt; ++ i) {
                cout << iter.instances[i].ins_id << '\t' << iter.instances[i].GetProb() << endl;
            }
        }
    }
};

int main(int argc, char const *argv[]) {
    // Dataset dataset();
    // dataset.gen
#ifdef _LINUX_
    struct timeval start, end;
    gettimeofday(&start, NULL);
#endif
    map<int, double> results;
    Dataset D(10, 4);
    double l[2] = {1, 2};
    double r[2] = {1, 2};
    HyperBox R(2, l, r);
    D.LoadData(argv[1]);
    // cout << "finish data load.\n";
    // D.PrintData();
    D.LoopPreprocessing();
    // cout << "finish preprocessing.\n";
    D.LoopAlg(R, results);
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


