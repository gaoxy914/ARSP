#include "object.hpp"

class Object {
public:
    int obj_id;
    int cnt;
    InstanceBase* instances;
};

class Dataset {
    int m, n;
    vector<Object> objects;

public:
    Dataset(const int& m) {
        n = 0;
        this->m = m;
        objects.resize(m);
    }

    void LoadData(const char* path) {
        ifstream file((string(path) + to_string(m) + string("/cnt.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) { file >> objects[i].cnt; n += objects[i].cnt; }
        file.close();
        file.open((string(path) + to_string(m) + string("/instances.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) {
            objects[i].instances = new InstanceBase[objects[i].cnt];
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

    void TransAlg(const HyperBox& R, map<int, double> results) {
        vector<double*> extremes;
        R.GetExtremes(extremes);
        int dim = extremes.size();
        vector<InstanceBase> new_instances;
        new_instances.reserve(n);
        for (auto obj : objects) {
            for (int i = 0; i < obj.cnt; ++ i) {
                InstanceBase cur_ins = obj.instances[i], new_ins;
                new_ins.obj_id = cur_ins.obj_id;
                new_ins.ins_id = cur_ins.ins_id;
                new_ins.prob = cur_ins.prob;
                new_ins.coord = new double[dim];
                for (int j = 0; j < dim; ++ j) new_ins.coord[j] = dot(cur_ins.coord, extremes[j]);
                new_instances.push_back(new_ins);
            }
        }

        // build kd-tree

        // calculate skyline probability
           
    }
};