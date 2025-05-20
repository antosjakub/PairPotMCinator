#include "objects.h"
using namespace std;

namespace consts {
    gsl_vector_float * zero = gsl_vector_float_calloc(3);
    gsl_vector_float * e1 = gsl_vector_float_calloc(3);
    gsl_vector_float * e2 = gsl_vector_float_calloc(3);
    gsl_vector_float * e3 = gsl_vector_float_calloc(3);
    void init() {
        gsl_vector_float_set_basis(e1, 0);
        gsl_vector_float_set_basis(e2, 1);
        gsl_vector_float_set_basis(e3, 2);
    }
}

//void print_gsl_vector(gsl_vector_float * v) {
//    for (int i=0; i<v->size; ++i) {
//        printf("%f ", gsl_vector_float_get(v,i));
//    }
//    printf("\n");
//    return;
//}
template <typename T> void print_vector_inline(vector<T> v) {
    for (int i=0; i<v.size()-1; ++i) {
        cout << v.at(i) << ", ";
    }
    cout << v.at(v.size()-1) << "\n";
    return;
}
template <typename T> void print_vector(vector<T> v) {
    for (int i=0; i<v.size(); ++i) {
        cout << v.at(i) << endl;
    }
    return;
}
void print_list(list<string> * list) {
    for(string el : *list) {
        cout << el << endl;
    }
    return;
}

void merge_two_lists(list<string> * dest, list<string> * src) {
    dest->merge(*src);
    dest->unique();
    delete src;
    return;
}

float sum_vector_elements(vector<float> v) {
    float sum = 0;
    vector<float>::iterator iter;
    for (iter=v.begin(); iter!=v.end(); ++iter) {
        sum += *iter;
    }
    return sum;
}
// TD: include inside mol_part definition
void vector_substract_one(vector<float> &v) {
    vector<float>::iterator iter;
    for (iter=v.begin(); iter!=v.end(); ++iter) {
        *iter = *iter - 1;
    }
}



float arccos(float x) {
    return x > 1.0 ? 1.0 : acos(x);
}

gsl_vector_float * init_gsl_vector(vector<float> v_in) {
    gsl_vector_float * v_out = gsl_vector_float_calloc(v_in.size());
    for (int i=0; i<v_in.size(); ++i) {
        gsl_vector_float_set(v_out, i, v_in.at(i));
    }
    return v_out;
};
// creates a mol list with 2 different charities
list<Molecule*> prepare_mol_list(Molecule * mol, int n_mols_chirality_0, int n_mols_chirality_1) {
    Molecule * mol_cp;
    list<Molecule*> mol_list;
    for (int i=0; i<n_mols_chirality_0; ++i) {
        mol_cp = new Molecule(mol->n_atoms);
        *mol_cp = *mol;
        mol_list.push_back(mol_cp);
    };
    for (int i=0; i<n_mols_chirality_1; ++i) {
        mol_cp = new Molecule(mol->n_atoms);
        *mol_cp = *mol;
        reverse_mol(mol_cp);
        mol_list.push_back(mol_cp);
    };
    return mol_list;
}
list<Molecule*> create_mol_list_from_mol(Molecule * mol, int n_mols) {
    Molecule * mol_cp;
    list<Molecule*> mol_list;
    for (int i=0; i<n_mols; ++i) {
        mol_cp = new Molecule(mol->n_atoms);
        *mol_cp = *mol;
        mol_list.push_back(mol_cp);
    };
    return mol_list;
}
list<Molecule*> create_mol_list_from_list(list<Molecule*> original) {
    Molecule * mol_cp;
    list<Molecule*> mol_list_cp;
    list<Molecule*>::iterator iter;
    for (iter=original.begin(); iter!=original.end(); ++iter) {
        mol_cp = new Molecule((*iter)->n_atoms);
        *mol_cp = *(*iter);
        mol_list_cp.push_back(mol_cp);
    };
    return mol_list_cp;
}
void cp_to_mol_list_from_list(list<Molecule*> destination, list<Molecule*> source) {
    list<Molecule*>::iterator iter_dest = destination.begin();
    list<Molecule*>::iterator iter_src = source.begin();
    for (iter_dest; iter_dest!=destination.end(); ++iter_dest, ++iter_src) {
        *(*iter_dest) = *(*iter_src);
    }
}
void append_mols_to_mol_list(list<Molecule*> & mol_list, Molecule * mol, int n_mols) {
    Molecule * mol_cp;
    for (int i=0; i<n_mols; ++i) {
        mol_cp = new Molecule(mol->n_atoms);
        *mol_cp = *mol;
        mol_list.push_back(mol_cp);
    };
}



int return_element_index(list<string> * keys, string key) {
    int i=0;
    int el_found=false;
    list<string>::iterator iter;
    for (iter=keys->begin(); iter!=keys->end(); ++iter) {
        if (*iter == key) {
            el_found=true;
            break;
        } else {
            ++i;
        }
    }
    if (el_found) {
        return i;
    } else {
        // throw error
        return -1;
    }
}
void populate_el_index(AtomHolder * holder, list<string> * used_elements) {
    for (int i=0; i<holder->n_atoms; ++i) {
        (holder->atom+i)->el_index = return_element_index(used_elements, (holder->atom+i)->element);
    }
}


Surface * return_cell_multiple(Surface * surf, int v1_coeff, int v2_coeff) {

    Surface * cell = new Surface(v1_coeff*v2_coeff*surf->n_atoms);
    cell->v1 = gsl_vector_float_calloc(3);
    cell->v2 = gsl_vector_float_calloc(3);

    Atom * surf_atom;
    gsl_vector_float * w1 = gsl_vector_float_calloc(3);
    gsl_vector_float * w2 = gsl_vector_float_calloc(3);
    gsl_vector_float * w = gsl_vector_float_calloc(3);
    gsl_vector_float * w_atom = gsl_vector_float_calloc(3);
    float x,y,z;
    int k=0;
    for (int i=0; i<v1_coeff; ++i) {
        gsl_vector_float_set_zero(w2);
        for (int j=0; j<v2_coeff; ++j) {
            gsl_vector_float_add(w,w1);
            gsl_vector_float_add(w,w2);
            for (int l=0; l<surf->n_atoms; ++l) {
                surf_atom = surf->atom+l;
                gsl_vector_float_memcpy(w_atom, w);
                gsl_vector_float_add(w_atom, surf_atom->pos);
                x = gsl_vector_float_get(w_atom, 0);
                y = gsl_vector_float_get(w_atom, 1);
                z = gsl_vector_float_get(w_atom, 2);
                (cell->atom+k)->fill(surf_atom->element, x,y,z, surf_atom->charge, surf_atom->hb_type, surf_atom->hb_acceptor);
                ++k;
            }
            gsl_vector_float_set_zero(w);
            gsl_vector_float_add(w2, surf->v2);
        }
        gsl_vector_float_add(w1, surf->v1);
    }
    gsl_vector_float_memcpy(cell->v1, surf->v1);
    gsl_vector_float_memcpy(cell->v2, surf->v2);
    gsl_vector_float_scale(cell->v1, v1_coeff);
    gsl_vector_float_scale(cell->v2, v2_coeff);

    return cell;
}



void substract_one_from_vector(vector<int> & v) {
    for (auto & el : v) {
        el = el - 1;
    }
    return;
}