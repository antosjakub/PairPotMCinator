#include "buckingham_pot_data.h"
#include <map>
#include <filesystem>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
using namespace std;


namespace Options {
    list<string> * used_elements;

    void init_surf(Surface * surf) {
        used_elements = new list<string>;
        merge_two_lists(Options::used_elements, surf->return_element_list());
        populate_el_index(surf, Options::used_elements);
    }
    void init_mol(Molecule * mol) {
        used_elements = new list<string>;
        merge_two_lists(Options::used_elements, mol->return_element_list());
        populate_el_index(mol, Options::used_elements);
    }
    void init_mol_surf(Molecule * mol, Surface * surf) {
        used_elements = new list<string>;
        merge_two_lists(Options::used_elements, mol->return_element_list());
        merge_two_lists(Options::used_elements, surf->return_element_list());
        populate_el_index(mol, Options::used_elements);
        populate_el_index(surf, Options::used_elements);
    }
    void init_mol_list_surf(list<Molecule*> mol_list, Surface * surf) {
        used_elements = new list<string>;

        list<Molecule*>::iterator iter;
        for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
            merge_two_lists(Options::used_elements, (*iter)->return_element_list());
        }
        merge_two_lists(Options::used_elements, surf->return_element_list());

        for (iter=mol_list.begin(); iter!=mol_list.end(); ++iter) {
            populate_el_index(*iter, Options::used_elements);
        }
        populate_el_index(surf, Options::used_elements);
    }
};



namespace PairPotentials {


class Electrostatic {
    public:
    float const COULOMB_CONSTANT = 14.3996;
    float return_en(Atom * atom_1, Atom * atom_2, float r) {
        return COULOMB_CONSTANT * atom_1->charge * atom_2->charge / r;
    }
};


class Buckingham {
    gsl_matrix_float ** buckingham_params;
    // data used in pauli & london potentials (together called a buckingham potential)
    public:
    Buckingham() {
        buckingham_params = new gsl_matrix_float*[3];
        buckingham_params[0] = init_buckingham_matrix(0);
        buckingham_params[1] = init_buckingham_matrix(1);
        buckingham_params[2] = init_buckingham_matrix(2);
    };
    private:
    list<string> * return_buckingham_elements() {
        list<string> * list_elements = new list<string>;
        map<string,vector<float>>::iterator iter;
        for (iter=buckingham_pot_data.begin(); iter!=buckingham_pot_data.end(); ++iter) {
            list_elements->push_back(iter->first);
        };
        return list_elements;
    }
    gsl_matrix_float * init_buckingham_matrix(int layer) {
        int n = (*Options::used_elements).size();
        // find values and copy them to vector
        vector<float> vals(n,0);
        map<string,vector<float>>::iterator iter_map;
        int i = -1;
        for(string key : *Options::used_elements) {
            iter_map = (buckingham_pot_data).find(key);
            vals.at(++i) = iter_map->second.at(layer);
        }
        // create matrix from vector values
        gsl_matrix_float * A = gsl_matrix_float_calloc(n,n);
        float a;
        for (int i=0; i<n; ++i) {
            for (int j=0; j<n; ++j) {
                a = ( i==j ? vals.at(i) : sqrt(vals.at(i)*vals.at(j)) );
                gsl_matrix_float_set(A, i,j, a);
            }
        }
        return A;
    };
    public:
    float get_param(int i, int j, int layer) {
        return gsl_matrix_float_get(buckingham_params[layer], i,j);
    }
    float return_en_pauli(Atom * atom_1, Atom * atom_2, float r) {
        // a*exp(-b*r)
        return get_param(atom_1->el_index,atom_2->el_index,0) * exp(-1*get_param(atom_1->el_index,atom_2->el_index,1)*r);
    }
    float return_en_london(Atom * atom_1, Atom * atom_2, float r) {
        // -c/r^6
        return r > 1.5 ? -1*get_param(atom_1->el_index,atom_2->el_index,2)/pow(r,6) : 0.0;
    }
    float return_en(Atom * atom_1, Atom * atom_2, float r) {
        return return_en_pauli(atom_1, atom_2, r) + return_en_london(atom_1, atom_2, r);
    }
};


class HydrogenBonding {
    public:
    float d=0.3, a=1.631, R=2.347;
    private:
    float return_en_hb_1d(Atom * atom_1, Atom * atom_2, float r, AtomHolder * mol_1, gsl_vector_float * offset = consts::zero, float offset_mult = 0.0) {
        gsl_vector_float * v_21 = gsl_vector_float_calloc(3);
        gsl_vector_float * v_31 = gsl_vector_float_calloc(3);
        float en;
        if (atom_1->hb_type == "H" && atom_2->hb_type == "D") {
            float theta;
            Atom * atom_3 = (mol_1->atom + atom_1->hb_acceptor);
            gsl_vector_float_memcpy(v_21, atom_2->pos);
            gsl_vector_float_memcpy(v_31, atom_3->pos);
            gsl_vector_float_sub(v_21, atom_1->pos);
            gsl_vector_float_sub(v_31, atom_1->pos);
            //gsl_vector_float_add(v_21, offset);
            gsl_blas_saxpy(offset_mult, offset, v_21);
            gsl_blas_sdot(v_21, v_31, &theta);
            theta /= gsl_blas_snrm2(v_21) * gsl_blas_snrm2(v_31);
            theta = arccos(theta);
            en = d * (exp(-2*a*(r-R))-2*exp(-a*(r-R))) * exp(-pow(theta/M_PI-1, 2));
        } else {
            en = 0;
        }
        gsl_vector_float_free(v_21);
        gsl_vector_float_free(v_31);
        return en;
    }
    public:
    float return_en(Atom * atom_1, Atom * atom_2, float r, AtomHolder * mol_1, AtomHolder * mol_2, gsl_vector_float * offset = consts::zero) {
        return return_en_hb_1d(atom_1, atom_2, r, mol_1, offset, 1.0) + 
               return_en_hb_1d(atom_2, atom_1, r, mol_2, offset, -1.0);
    }
};



PairPotentials::Electrostatic * elstat;
PairPotentials::Buckingham * buckingham;
PairPotentials::HydrogenBonding * hb;

void init() {
    elstat = new PairPotentials::Electrostatic;
    buckingham = new PairPotentials::Buckingham;
    hb = new PairPotentials::HydrogenBonding;
}

};


// make as an namespace -> globally available
void calc_en_mol_within(Molecule * mol, gsl_vector_float * en_vec) {
    Atom * atom_1;
    Atom * atom_2;
    gsl_vector_float * r_vec = gsl_vector_float_alloc(3);
    float r = 0;
    for (int i=0; i<mol->n_atoms; ++i) {
        atom_1 = mol->atom+i;
        for (int j=i+1; j<mol->n_atoms; ++j) {
            atom_2 = mol->atom+j;
            gsl_vector_float_memcpy(r_vec, atom_2->pos);
            gsl_vector_float_sub(r_vec, atom_1->pos);
            r = gsl_blas_snrm2(r_vec);
            gsl_vector_float_set(en_vec, 1, gsl_vector_float_get(en_vec,1) + PairPotentials::elstat->return_en(atom_1, atom_2, r));
            gsl_vector_float_set(en_vec, 2, gsl_vector_float_get(en_vec,2) + PairPotentials::buckingham->return_en_pauli(atom_1, atom_2, r));
            gsl_vector_float_set(en_vec, 3, gsl_vector_float_get(en_vec,3) + PairPotentials::buckingham->return_en_london(atom_1, atom_2, r));
            gsl_vector_float_set(en_vec, 0, gsl_vector_float_get(en_vec,0) + PairPotentials::hb->return_en(atom_1, atom_2, r, mol, mol));
        }
    }
    gsl_vector_float_free(r_vec);
    return;
}
void calc_en_mol_surf(Molecule * mol, Surface * surf, gsl_vector_float * en_vec, gsl_vector_float * offset = consts::zero) {
    Atom * atom_1;
    Atom * atom_2;
    gsl_vector_float * r_vec = gsl_vector_float_alloc(3);
    float r = 0;
    for (int i=0; i<mol->n_atoms; ++i) {
        atom_1 = mol->atom+i;
        for (int j=0; j<surf->n_atoms; ++j) {
            atom_2 = surf->atom+j;
            gsl_vector_float_memcpy(r_vec, atom_2->pos);
            gsl_vector_float_sub(r_vec, atom_1->pos);
            gsl_vector_float_add(r_vec, offset);
            r = gsl_blas_snrm2(r_vec);
            gsl_vector_float_set(en_vec, 2, gsl_vector_float_get(en_vec,2) + PairPotentials::buckingham->return_en_pauli(atom_1, atom_2, r));
            gsl_vector_float_set(en_vec, 3, gsl_vector_float_get(en_vec,3) + PairPotentials::buckingham->return_en_london(atom_1, atom_2, r));
        }
    }
    gsl_vector_float_free(r_vec);
    return;
}
void calc_en_mol_mol(Molecule * mol_1, Molecule * mol_2, gsl_vector_float * en_vec, gsl_vector_float * offset = consts::zero) {
    Atom * atom_1;
    Atom * atom_2;
    gsl_vector_float * r_vec = gsl_vector_float_alloc(3);
    float r = 0;
    for (int i=0; i<mol_1->n_atoms; ++i) {
        atom_1 = mol_1->atom+i;
        for (int j=0; j<mol_2->n_atoms; ++j) {
            atom_2 = mol_2->atom+j;
            gsl_vector_float_memcpy(r_vec, atom_2->pos);
            gsl_vector_float_sub(r_vec, atom_1->pos);
            gsl_vector_float_add(r_vec, offset);
            r = gsl_blas_snrm2(r_vec);
            gsl_vector_float_set(en_vec, 1, gsl_vector_float_get(en_vec,1) + PairPotentials::elstat->return_en(atom_1, atom_2, r));
            gsl_vector_float_set(en_vec, 2, gsl_vector_float_get(en_vec,2) + PairPotentials::buckingham->return_en_pauli(atom_1, atom_2, r));
            gsl_vector_float_set(en_vec, 3, gsl_vector_float_get(en_vec,3) + PairPotentials::buckingham->return_en_london(atom_1, atom_2, r));
            gsl_vector_float_set(en_vec, 0, gsl_vector_float_get(en_vec,0) + PairPotentials::hb->return_en(atom_1, atom_2, r, mol_1, mol_2, offset));
        }
    }
    gsl_vector_float_free(r_vec);
    return;
}
void calc_en_mol_list(list <Molecule*> mol_list, gsl_vector_float * en_vec) {
    list<Molecule*>::iterator iter_1;
    list<Molecule*>::iterator iter_2;
    Atom * atom_1;
    Atom * atom_2;
    gsl_vector_float * r_vec = gsl_vector_float_alloc(3);
    float r = 0;
    for (iter_1=mol_list.begin(); iter_1!=mol_list.end(); ++iter_1) {
        iter_2=iter_1;
        for (++iter_2; iter_2!=mol_list.end(); ++iter_2) {
            calc_en_mol_mol(*iter_1, *iter_2, en_vec);
        }
    }
    gsl_vector_float_free(r_vec);
    return;
}


void descent_mol(Molecule * mol, Surface * surface) {

    float initial_height = 6.0;
    float step_size = -0.1;
    gsl_vector_float * v = gsl_vector_float_calloc(3);
    gsl_vector_float_set(v, 2, initial_height);
    mol->move(v);
    //write_xyz_single_surface(mol, surface, "../xyz_files/playground.xyz", "w");

    float en_new = 0, en_old = 0;
    gsl_vector_float * en_vec = gsl_vector_float_calloc(4);
    calc_en_mol_surf(mol, surface, en_vec);
    en_old = gsl_vector_float_sum(en_vec);
    gsl_vector_float_set(v, 2, step_size);
    int n_steps = 0;
    //std::cout << n_steps << ": " << en_old << endl;
    while (true) {
        mol->move(v);
        gsl_vector_float_set_zero(en_vec);
        calc_en_mol_surf(mol, surface, en_vec);
        en_new = gsl_vector_float_sum(en_vec);
        if (en_new < en_old) {
            en_old = en_new;
            ++n_steps;
            //std::cout << n_steps << ": " << en_old << endl;
            //write_xyz_single_surface(mol, surface, "../xyz_files/playground.xyz", "a");
        } else {
            //std::cout << -1 << ": " << en_new << endl;
            std::cout << "descent_mol ended at height = " << initial_height + step_size*n_steps << ", with energy = " << en_old << endl;
            break;
        }
    }
    return;
}