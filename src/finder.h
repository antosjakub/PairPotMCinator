#include "energy.h"
#include <thread>
#include <limits>


gsl_matrix_float * compute_basis_inverse(gsl_vector_float * v1, gsl_vector_float * v2) {
    float det_A = gsl_vector_float_get(v1, 0) * gsl_vector_float_get(v2, 1)
                - gsl_vector_float_get(v2, 0) * gsl_vector_float_get(v1, 1);
    gsl_matrix_float * A_inv = gsl_matrix_float_calloc(3,3);
    gsl_matrix_float_set(A_inv, 0,0, gsl_vector_float_get(v2, 1));
    gsl_matrix_float_set(A_inv, 0,1, -1*gsl_vector_float_get(v2, 0));
    gsl_matrix_float_set(A_inv, 1,0, -1*gsl_vector_float_get(v1, 1));
    gsl_matrix_float_set(A_inv, 1,1, gsl_vector_float_get(v1, 0));
    gsl_matrix_float_scale(A_inv, 1.0/det_A);
    gsl_matrix_float_set(A_inv, 2,2, 1);
    return A_inv;
};

map<string, list<string>> log_source = {
    {"txt", {
        "cycle_min",
        "stages_min",
        "all_min",
    }},
    {"xyz", {
        "cycle_min_cell",
        "cycle_min_periodic",
        "stages_min_cell",
        "stages_min_periodic",
        "all_min_cell",
        "all_min_periodic"
    }}
};
class LogObject {
    public:
    bool use;
    string path;
    ofstream file;
    LogObject() : use(false), path("") {}
    LogObject(bool log_use, string log_path) {
        use = log_use;
        path = log_path;
    }
    ~LogObject() {
        if (file.is_open()) {
            file.close();
        }
    }
};
//class LoggingObject {
//    public:
//    map<string, map<string, LogObject>> logging_metainfo;
//    LoggingObject(map<string, map<string, LogObject>> logging_metainfo_in) {
//        logging_metainfo = logging_metainfo_in;
//    }
//    void new_cycle(int index) {
//        std::cout << "CYCLE: " << index << std::endl;
//        logging_metainfo["txt"]["cycle_min"].file << "CYCLE: " << index << std::endl;
//        logging_metainfo["txt"]["stages_min"].file << "CYCLE: " << index << std::endl;
//        logging_metainfo["txt"]["all_min"].file << "CYCLE: " << index << std::endl;
//    }
//    void end_cycle() {
//        std::cout << "END OF CYCLE" << std::endl;
//        //logging_metainfo["xyz"]["cycle_min"].file;
//    }
//};
//class LogingObject {
//    public:
//    map<string, map<string, LogObject>> logging_metainfo;
//    LogingObject(map<string, map<string, LogObject>> logging_metainfo_in) {
//        logging_metainfo = logging_metainfo_in;
//    }
//    void new_(string log_lvl, int index) {
//        std::cout << log_lvl << ": " << index << std::endl;
//        str[3] = ["cycle", "stage", "loop"];
//        //"cycle", 0
//        //"stage", 1
//        //"loop", 2
//        for(int i=0; i<=log_lvl, i++) {
//            logging_metainfo["txt"][str[i]].file << "CYCLE: " << index << std::endl;
//        }
//        logging_metainfo["txt"]["stages_min"].file << "CYCLE: " << index << std::endl;
//        logging_metainfo["txt"]["all_min"].file << "CYCLE: " << index << std::endl;
//    }
//    void end_() {
//        std::cout << "END OF STAGE" << std::endl;
//        //logging_metainfo["xyz"]["cycle_min"].file;
//    }
//};


// use in every stage for each group
// use DOF: T=int, multipliers: T=float, bind: T=int

class Stage {
    public:
    string name;
    int n_loops;
    int n_trials;
    map<string, map<string, int>> bind;
    map<string, map<string, float>> params;
    public:
    Stage(string name_in, int n_loops_in, int n_trials_in) {
        name = name_in;
        n_loops = n_loops_in;
        n_trials = n_trials_in;
    }
    void print_out() {
        cout << name << endl;
        cout << "n_loops: " << n_loops << endl;
        cout << "n_trials: " << n_trials << endl;
    }
};


class Finder {
    public:

    // constructor
    Surface * cell;
    int n_cycles;

    map<string, map<string, LogObject*>> log;

    // preallocate useful variables
    float v1_norm, v2_norm;
    gsl_vector_float * v1, * v2;
    gsl_vector_float * v1_unit, * v2_unit;
    gsl_vector_float * en_vec_within_original;
    float en_min;
    float en_found;

    Molecule * mol;
    list<Molecule*> mol_list_original;
    list<Molecule*> mol_list_curr;
    list<Molecule*> mol_list_cp;
    list<Molecule*> mol_list_winner;

    //private:
    gsl_matrix_float *A_inv;
    vector<int> mults;

    map<float, list<Molecule*>> thread_results;
    map<float, list<Molecule*>>::iterator thread_iter;
    map<string,int> mol_group_counts;

    int seed_value;


    //Finder(Surface * cell_in, list<Molecule*> & mol_list_in, list<MoleculeGroup*> & mol_group_list_in, int n_cycles_in, int seed) {
    Finder(Surface * cell_in, list<Molecule*> & mol_list_in, map<string, int> group_counts, int n_cycles_in, int seed, map<string, map<string, LogObject*>> log_in) {
        log = log_in;
        mol_group_counts = group_counts; // info about passed molecule groups

        cell = cell_in;
        mol_list_original = mol_list_in;
        n_cycles = n_cycles_in;

        // calc en within each mol and sum it to a vector
        en_vec_within_original = gsl_vector_float_calloc(4);
        list<Molecule*>::iterator iter;
        for (iter=mol_list_original.begin(); iter!=mol_list_original.end(); ++iter) {
            calc_en_mol_within(*iter, en_vec_within_original);
        }

        mol_list_curr = create_mol_list_from_list(mol_list_original);
        mol_list_cp = create_mol_list_from_list(mol_list_original);
        mol_list_winner = create_mol_list_from_list(mol_list_original);

        //rand = gsl_rng_alloc(gsl_rng_taus);
        //gsl_rng_set(rand, seed);
        seed_value = seed;

        // setting up basis and normed basis
        v1 = cell->v1;
        v2 = cell->v2;
        v1_unit = gsl_vector_float_alloc(3);
        v2_unit = gsl_vector_float_alloc(3);
        gsl_vector_float_memcpy(v1_unit, v1);
        gsl_vector_float_memcpy(v2_unit, v2);
        v1_norm = gsl_blas_snrm2(v1);
        v2_norm = gsl_blas_snrm2(v2);
        gsl_vector_float_scale(v1_unit, 1/v1_norm);
        gsl_vector_float_scale(v2_unit, 1/v2_norm);

        mults = {1,0, 1,1, 0,1, -1,1, -1,0, -1,-1, 0,-1, 1,-1};
        A_inv = compute_basis_inverse(v1,v2);
    }

    void variate(list<Molecule*> mol_list, map<string, int> bind, map<string, float> params, gsl_rng * rand) {
        float phi, a1, a2, a3;
        gsl_vector_float * w = gsl_vector_float_calloc(3);
        gsl_vector_float * unit_vector = gsl_vector_float_calloc(3);
        list<Molecule*>::iterator iter;

        if (params["rotate_e3"] > 0) {
            if (bind["rotate_e3"] == 1) {
                phi = params["rotate_e3"] * M_PI*(2*gsl_rng_uniform(rand)-1);
                for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
                    (*iter)->rotate(consts::e3, phi);
                }
            } else {
                for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
                    phi = params["rotate_e3"] * M_PI*(2*gsl_rng_uniform(rand)-1);
                    (*iter)->rotate(consts::e3, phi);
                }
            }
        };
        if (params["rotate_3D"] > 0) {
            if (bind["rotate_3D"] == 1) {
                phi = params["rotate_3D"] * M_PI*(2*gsl_rng_uniform(rand)-1);
                set_random_axis(unit_vector, rand);
                for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
                    (*iter)->rotate(unit_vector, phi);
                }
            } else {
                for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
                    phi = params["rotate_3D"] * M_PI*(2*gsl_rng_uniform(rand)-1);
                    set_random_axis(unit_vector, rand);
                    (*iter)->rotate(unit_vector, phi);
                }
            }
        };
        if (params["translate_2D"] > 0) {
            if (bind["translate_2D"] == 1) {
                a1 = params["translate_2D"] * (2*gsl_rng_uniform(rand)-1);
                a2 = params["translate_2D"] * (2*gsl_rng_uniform(rand)-1);
                gsl_vector_float_set_zero(w);
                gsl_blas_saxpy(a1, v1_unit, w);
                gsl_blas_saxpy(a2, v2_unit, w);
                for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
                    (*iter)->move(w);
                }
            } else {
                for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
                    a1 = params["translate_2D"] * (2*gsl_rng_uniform(rand)-1);
                    a2 = params["translate_2D"] * (2*gsl_rng_uniform(rand)-1);
                    gsl_vector_float_set_zero(w);
                    gsl_blas_saxpy(a1, v1_unit, w);
                    gsl_blas_saxpy(a2, v2_unit, w);
                    (*iter)->move(w);
                }
            }
        };
        if (params["translate_e3"] > 0) {
            if (bind["translate_e3"] == 1) {
                a3 = params["translate_e3"] * (2*gsl_rng_uniform(rand)-1);
                gsl_vector_float_set_zero(w);
                gsl_vector_float_set(w, 2, a3);
                for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
                    (*iter)->move(w);
                }
            } else {
                for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
                    a3 = params["translate_e3"] * (2*gsl_rng_uniform(rand)-1);
                    gsl_vector_float_set_zero(w);
                    gsl_vector_float_set(w, 2, a3);
                    (*iter)->move(w);
                }
            }
        };
        if (params["deform"] > 0) {
            if (bind["deform"] == 1) {
                int n_aa_deforms = ((*mol_list.begin())->parts_around_axis).size();
                vector<float> rands_aa_deform(n_aa_deforms);
                for (int i=0; i<n_aa_deforms; ++i) {
                    phi = params["deform"] * M_PI*(2*gsl_rng_uniform(rand)-1);
                    //phi = M_PI*gsl_rng_uniform_int(rand, 2); // rand int gives {0,1}
                    rands_aa_deform.at(i) = phi;
                }

                int n_fe_deforms = ((*mol_list.begin())->parts_free_end).size();
                vector<float> rands_fe_deform(n_fe_deforms);
                vector<gsl_vector_float*> rands_fe_deform_vec(n_fe_deforms);
                for (int i=0; i<n_fe_deforms; ++i) {
                    phi = params["deform"] * M_PI*(2*gsl_rng_uniform(rand)-1);
                    rands_fe_deform.at(i) = phi;
                    // td
                    gsl_vector_float * rand_vec = gsl_vector_float_calloc(3);
                    set_random_axis(rand_vec, rand);
                    rands_fe_deform_vec.at(i) = rand_vec;
                }

                int i;
                for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
                    Atom * atom = (*iter)->atom;
                    i = 0; for (auto part : (*iter)->parts_around_axis) {
                        part.rotate(rands_aa_deform.at(i), atom);
                        ++i;
                    }
                    i = 0; for (auto part : (*iter)->parts_free_end) {
                        part.rotate(rands_fe_deform.at(i), atom, rands_fe_deform_vec.at(i));
                        ++i;
                    }
                }
                for (int i=0; i<n_fe_deforms; ++i) {
                    gsl_vector_float_free(rands_fe_deform_vec.at(i));
                }

            } else {
                for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
                    Atom * atom = (*iter)->atom;
                    for (auto part : (*iter)->parts_around_axis) {
                        phi = params["deform"] * M_PI*(2*gsl_rng_uniform(rand)-1);
                        part.rotate(phi, atom);
                    }
                    for (auto part : (*iter)->parts_free_end) {
                        phi = params["deform"] * M_PI*(2*gsl_rng_uniform(rand)-1);
                        set_random_axis(unit_vector, rand);
                        part.rotate(phi, atom, unit_vector);
                    }
                }
            }
        }
        for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
            periodically_cut(*iter);
        }
        gsl_vector_float_free(w);
        gsl_vector_float_free(unit_vector);
    }

    void variate_wrapper(float en_min_global, list<Molecule*> mol_list, map<string, map<string, int>> bind, map<string, map<string, float>> params, gsl_rng * rand) {
        int start_i = 0;
        int end_i = 0;
        list<Molecule*>::iterator iter;
        list<Molecule*> mol_list_copy = create_mol_list_from_list(mol_list);
        for (auto group_counts : mol_group_counts) {
            string name = group_counts.first;
            int n_mols = group_counts.second;
            end_i = start_i + n_mols;
            auto start = next(mol_list_copy.begin(), start_i);
            auto end = next(mol_list_copy.begin(), end_i);
            start_i = end_i;
            list<Molecule*> mol_list_group(start, end);
            variate(mol_list_group, bind[name], params[name], rand);
        }

        // en
        float en = calc_en_mol_list_periodic(mol_list_copy);
        //cout << "---- en_found: " << en << ", en_min_global: " << en_min_global << endl;
        //write_xyz_list_surface(mol_list_copy, cell, "playground.xyz", "a", "craziness");
        if (en < en_min_global) {
            thread_results[en] = mol_list_copy;
        } else {
            mol_list_copy.clear();
        }
    }



    float variate_threads(int n_trials, float en_min_global, list<Molecule*> mol_list, map<string, map<string, int>> bind,  map<string, map<string, float>> params, vector<gsl_rng*> random_generators) {
        thread_results.clear();
        float en_min_local = en_min_global;
        vector<std::thread> threads(n_trials);
        for (int i=0; i<n_trials; ++i) {
            threads.at(i) = std::thread(&Finder::variate_wrapper, this, en_min_local, mol_list, bind, params, random_generators.at(i));
        }
        for (auto& th : threads) {
            th.join();
        }

        thread_iter = thread_results.begin();
        if (thread_iter != thread_results.end()) {
            en_found = (*thread_iter).first;
            if (en_found < en_min_local) {
                en_min_local = en_found;
                cp_to_mol_list_from_list(mol_list_winner, (*thread_iter).second);
            }
        }

        threads.clear();
        return en_min_local;
    }

    void prepare_for_new_cycle() {
        en_min = std::numeric_limits<int>::max();
        cp_to_mol_list_from_list(mol_list_curr, mol_list_original);
    }
    
    //vector<gsl_rng*> append_random_generators(int n_trials) {


    void run_scatter(Stage stage, string stage_string) {
        // vector<gsl_rng*> random_generators = append_random_generators(n_trials)
        vector<gsl_rng*> random_generators(stage.n_trials);
        for (int j=0; j<stage.n_trials; ++j) {
            random_generators.at(j) = gsl_rng_alloc(gsl_rng_taus);
            gsl_rng_set(random_generators.at(j), ++seed_value);
        }
        en_found = variate_threads(stage.n_trials, en_min, mol_list_curr, stage.bind, stage.params, random_generators);
        std::cout << stage_string << " | en: " << en_found << endl;
        if (en_found < en_min) {
            en_min = en_found;
            cp_to_mol_list_from_list(mol_list_curr, mol_list_winner);
            if (log["txt"]["all_min"]->use) log["txt"]["all_min"]->file << stage_string << " | en: " << en_min << endl;
            if (log["xyz"]["all_min_cell"]->use) write_xyz_list_surface(mol_list_curr, cell, log["xyz"]["all_min_cell"]->path, "a", stage_string+", scatter");
            if (log["xyz"]["all_min_periodic"]->use) write_xyz_list_surface_periodic(cell->v1, cell->v2, mol_list_curr, cell, log["xyz"]["all_min_periodic"]->path, "a", stage_string+", scatter");

        }
        random_generators.clear();
        return;
    }

    void run_stage(Stage stage, string stage_string) {
        // vector<gsl_rng*> random_generators = append_random_generators(n_trials)
        vector<gsl_rng*> random_generators(stage.n_trials);
        for (int j=0; j<stage.n_trials; ++j) {
            random_generators.at(j) = gsl_rng_alloc(gsl_rng_taus);
            gsl_rng_set(random_generators.at(j), ++seed_value);
        }
        for (int j=0; j<stage.n_loops; ++j) {
            en_found = variate_threads(stage.n_trials, en_min, mol_list_curr, stage.bind, stage.params, random_generators);
            std::cout << stage_string << ", loop: " << j+1 << " | en: " << en_found << endl;
            if (en_found < en_min) {
                en_min = en_found;
                cp_to_mol_list_from_list(mol_list_curr, mol_list_winner);
                if (log["txt"]["all_min"]->use) log["txt"]["all_min"]->file << stage_string+", loop: " << j+1 << " | en: " << en_min << endl;
                if (log["xyz"]["all_min_cell"]->use) write_xyz_list_surface(mol_list_curr, cell, log["xyz"]["all_min_cell"]->path, "a", stage_string+", loop: "+to_string(j+1));
                if (log["xyz"]["all_min_periodic"]->use) write_xyz_list_surface_periodic(cell->v1, cell->v2, mol_list_curr, cell, log["xyz"]["all_min_periodic"]->path, "a", stage_string+", loop: "+to_string(j+1));
            }
        }
        random_generators.clear();
        return;
    }


    float calc_en_mol_list_periodic(list<Molecule*> mol_list_local) {
        list<Molecule*>::iterator iter_1;
        list<Molecule*>::iterator iter_2;

        gsl_vector_float * en_vec_mol_surf = gsl_vector_float_calloc(4);
        gsl_vector_float * en_vec_same_mol = gsl_vector_float_calloc(4);
        gsl_vector_float * en_vec_diff_mol = gsl_vector_float_calloc(4);
        gsl_vector_float * en_vec_within = gsl_vector_float_calloc(4);
        gsl_vector_float * en_vec = gsl_vector_float_calloc(4);

        gsl_vector_float *offset = gsl_vector_float_calloc(3);

        // calc en in center cell
        calc_en_mol_list(mol_list_local, en_vec_diff_mol);
        for (iter_1=mol_list_local.begin(); iter_1!=mol_list_local.end(); ++iter_1) {
            calc_en_mol_surf(*iter_1, cell, en_vec_mol_surf, consts::zero);
            calc_en_mol_within(*iter_1, en_vec_within);
        }
        gsl_vector_float_sub(en_vec_within, en_vec_within_original);

        // calc en in other cells
        for (int i=0; i<(mults.size()/2); ++i) {
            gsl_vector_float_set_zero(offset);
            gsl_blas_saxpy(mults.at(2*i), v1, offset);
            gsl_blas_saxpy(mults.at(2*i+1), v2, offset);
            for (iter_1=mol_list_local.begin(); iter_1!=mol_list_local.end(); ++iter_1) {
                calc_en_mol_surf(*iter_1, cell, en_vec_mol_surf, offset);
                if (i<4) {calc_en_mol_mol(*iter_1, *iter_1, en_vec_same_mol, offset);}
                iter_2=iter_1;
                for (++iter_2; iter_2!=mol_list_local.end(); ++iter_2) {
                    calc_en_mol_mol(*iter_1, *iter_2, en_vec_diff_mol, offset);
                }
            }
        }
        gsl_vector_float_scale(en_vec_same_mol, 2);

        gsl_vector_float_add(en_vec, en_vec_same_mol);
        gsl_vector_float_add(en_vec, en_vec_diff_mol);
        gsl_vector_float_add(en_vec, en_vec_mol_surf);
        gsl_vector_float_add(en_vec, en_vec_within);
        //print_gsl_vector(en_vec_same_mol);
        //print_gsl_vector(en_vec_diff_mol);
        //print_gsl_vector(en_vec_mol_surf);
        //print_gsl_vector(en_vec);
        
        gsl_vector_float_free(en_vec_mol_surf);
        gsl_vector_float_free(en_vec_same_mol);
        gsl_vector_float_free(en_vec_diff_mol);
        gsl_vector_float_free(en_vec_within);

        float energy = gsl_vector_float_sum(en_vec);
        gsl_vector_float_free(en_vec);
        gsl_vector_float_free(offset);
        return energy;
    }


    void periodically_cut(Molecule * mol) {
        gsl_vector_float *v = gsl_vector_float_calloc(3);
        gsl_vector_float *u = gsl_vector_float_calloc(3);
        // take molecule, get the coefficients of the center of mass vector in the basis {v1, v2, e3}
        // if it lies ouside of the cell given by (-0.5,0.5)v1 x (-0.5,0.5)v2, then move it inside...
        gsl_blas_sgemv(CblasNoTrans, 1.0, A_inv, mol->pos, 0.0, u);
        if (gsl_vector_float_get(u,0) > 0.5) {
            gsl_vector_float_memcpy(v, v1);
            gsl_vector_float_scale(v, -1);
            mol->move(v);
        } else if (gsl_vector_float_get(u,0) < -0.5) {
            gsl_vector_float_memcpy(v, v1);
            mol->move(v);
        }
        if (gsl_vector_float_get(u,1) > 0.5) {
            gsl_vector_float_memcpy(v, v2);
            gsl_vector_float_scale(v, -1);
            mol->move(v);
        } else if (gsl_vector_float_get(u,1) < -0.5) {
            gsl_vector_float_memcpy(v, v2);
            mol->move(v);
        }
        gsl_vector_float_free(u);
        gsl_vector_float_free(v);
        return;
    }

};