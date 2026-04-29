#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <cmath>

using namespace std;


void print_gsl_vector(gsl_vector_float * v) {
    for (int i=0; i<v->size; ++i) {
        printf("%f ", gsl_vector_float_get(v,i));
    }
    printf("\n");
    return;
}

void set_random_axis(gsl_vector_float * v, gsl_rng * rand) {
    gsl_vector_float_set(v, 0, (2*gsl_rng_uniform(rand)-1));
    gsl_vector_float_set(v, 1, (2*gsl_rng_uniform(rand)-1));
    gsl_vector_float_set(v, 2, (2*gsl_rng_uniform(rand)-1));
    gsl_vector_float_scale(v, 1.0/gsl_blas_snrm2(v));
    return;
}

class Atom {
    public:
        string element;
        int el_index;
        gsl_vector_float * pos;
        float charge;
        string hb_type;
        int hb_acceptor;
    Atom() {
        pos = gsl_vector_float_alloc(3);
    }
    Atom(Atom &t) {
        //cout << "Atom copy const called" << endl;
        pos = gsl_vector_float_alloc(3);
        gsl_vector_float_memcpy(pos, t.pos);
        element = t.element;
        el_index = t.el_index;
        charge = t.charge;
        hb_type = t.hb_type;
        hb_acceptor = t.hb_acceptor;
    }
    Atom& operator=(const Atom &t) {
        //cout << "Atom assign op called" << endl;
        gsl_vector_float_memcpy(pos, t.pos);
        element = t.element;
        el_index = t.el_index;
        charge = t.charge;
        hb_type = t.hb_type;
        hb_acceptor = t.hb_acceptor;
        return *this;
    }
    void fill(string element_ext, float x, float y, float z, float charge_ext, string hb_type_ext, int hb_acceptor_ext) {
        element = element_ext;
        gsl_vector_float_set(pos, 0, x);
        gsl_vector_float_set(pos, 1, y);
        gsl_vector_float_set(pos, 2, z);
        charge = charge_ext;
        hb_type = hb_type_ext;
        hb_acceptor = hb_acceptor_ext;
    }
};


void set_quaternion_matrix(const gsl_vector_float * axis, const float angle, gsl_matrix_float * A) {
    gsl_matrix_float_set_zero(A);
    gsl_matrix_float_set(A,0,1, -gsl_vector_float_get(axis, 2));
    gsl_matrix_float_set(A,0,2,  gsl_vector_float_get(axis, 1));
    gsl_matrix_float_set(A,1,0,  gsl_vector_float_get(axis, 2));
    gsl_matrix_float_set(A,1,2, -gsl_vector_float_get(axis, 0));
    gsl_matrix_float_set(A,2,0, -gsl_vector_float_get(axis, 1));
    gsl_matrix_float_set(A,2,1,  gsl_vector_float_get(axis, 0));
    gsl_matrix_float_scale(A, sin(angle));
    gsl_matrix_float_set(A,0,0, cos(angle));
    gsl_matrix_float_set(A,1,1, cos(angle));
    gsl_matrix_float_set(A,2,2, cos(angle));
    // A + a xy^T
    gsl_blas_sger(1-cos(angle), axis, axis, A); 
}
gsl_matrix_float * return_quaternion_matrix(const gsl_vector_float * axis, const float angle) {
    gsl_matrix_float * A = gsl_matrix_float_alloc(3,3);
    set_quaternion_matrix(axis, angle, A);
    return A;
}
void rotate_vec(gsl_vector_float * vec, const gsl_matrix_float * A, const gsl_vector_float * pos_0, gsl_vector_float * y) {
    // y = a Ax + b y
    gsl_vector_float_sub(vec, pos_0);
    gsl_blas_sgemv(CblasNoTrans, 1.0, A, vec, 0.0, y);
    gsl_vector_float_memcpy(vec, y);
    gsl_vector_float_add(vec, pos_0);
    return;
}


class AtomHolder {
    public:
        int n_atoms;
        Atom * atom;
        gsl_vector_float * pos;
    private:
        gsl_vector_float * __w;
        gsl_matrix_float * __A;
    public:
        explicit AtomHolder(int n_atoms_ext) {
            n_atoms = n_atoms_ext;
            atom = new Atom[n_atoms];
            pos = gsl_vector_float_calloc(3);
            __w = gsl_vector_float_calloc(3);
            __A = gsl_matrix_float_calloc(3,3);
        }
        AtomHolder(AtomHolder &t) {
            //cout << "copy constructor of AtomHolder called" << endl;
            n_atoms = t.n_atoms;
            atom = new Atom[n_atoms];
            for (int i=0; i<n_atoms; ++i) {
                *(atom+i) = *(t.atom+i);
            }
            pos = gsl_vector_float_alloc(3);
            __w = gsl_vector_float_alloc(3);
            __A = gsl_matrix_float_alloc(3,3);
            gsl_vector_float_memcpy(pos, t.pos);
            gsl_vector_float_memcpy(__w, t.__w);
            gsl_matrix_float_memcpy(__A, t.__A);
        };
        AtomHolder& operator=(const AtomHolder &t) {
            //cout << "assignment operator of AtomHolder called" << endl;
            n_atoms = t.n_atoms;
            for (int i=0; i<n_atoms; ++i) {
                *(atom+i) = *(t.atom+i);
            }
            gsl_vector_float_memcpy(pos, t.pos);
            gsl_vector_float_memcpy(__w, t.__w);
            gsl_matrix_float_memcpy(__A, t.__A);
            return *this;
        }
        list<string> * return_element_list();
        void set_center_of_mass(gsl_vector_float * v);
        void move(const gsl_vector_float * v);
        void set_pos(const gsl_vector_float * v);
        void center();
        void rotate(const gsl_vector_float * axis, const float angle);
};
list<string> * AtomHolder::return_element_list() {
    list<string> * element_list = new list<string>;
    string el;
    for (int i=0; i<n_atoms; ++i) {
        element_list->push_back((atom+i)->element);
    }
    element_list->sort();
    element_list->unique();
    return element_list;
}
void AtomHolder::set_center_of_mass(gsl_vector_float * v) {
    gsl_vector_float_set_zero(v);
    for (int i=0; i<n_atoms; ++i) {
        gsl_vector_float_add(v, (atom+i)->pos);
    }
    gsl_vector_float_scale(v, 1.0/n_atoms);
    return;
}
void AtomHolder::center() {
    set_center_of_mass(__w);
    for (int i=0; i<n_atoms; ++i) {
        gsl_vector_float_sub((atom+i)->pos, __w);
    }
    return;
}
void AtomHolder::rotate(const gsl_vector_float * axis, const float angle) {
    set_quaternion_matrix(axis, angle, __A);
    for (int i=0; i<n_atoms; ++i) {
        rotate_vec((atom+i)->pos, __A, pos, __w);
    };
    return;
}
void AtomHolder::move(const gsl_vector_float * v) {
    gsl_vector_float_add(pos, v);
    for (int i=0; i<n_atoms; ++i) {
        gsl_vector_float_add((atom+i)->pos, v);
    }
    return;
};
void AtomHolder::set_pos(const gsl_vector_float * v) {
    gsl_vector_float_memcpy(__w, v);
    gsl_vector_float_sub(__w, pos);
    for (int i=0; i<n_atoms; ++i) {
        gsl_vector_float_add((atom+i)->pos, __w);
    }
    return;
};


class MoleculePart {
    public:
        vector<int> atoms_indices;
        vector<int> axis_indices;
        gsl_vector_float * __rot_axis;
    private:
        gsl_vector_float * __w;
        gsl_matrix_float * __A;
    public:
    explicit MoleculePart(vector<int> axis_indices_inp, vector<int> atoms_indices_inp) {
        axis_indices = axis_indices_inp;
        atoms_indices = atoms_indices_inp;
        __rot_axis = gsl_vector_float_calloc(3);
        __w = gsl_vector_float_calloc(3);
        __A = gsl_matrix_float_calloc(3,3);
    }
    MoleculePart(const MoleculePart& other) {
        atoms_indices = other.atoms_indices;
        axis_indices = other.axis_indices;
        __rot_axis = gsl_vector_float_alloc(3);
        __w = gsl_vector_float_alloc(3);
        __A = gsl_matrix_float_alloc(3, 3);
        gsl_vector_float_memcpy(__rot_axis, other.__rot_axis);
        gsl_vector_float_memcpy(__w, other.__w);
        gsl_matrix_float_memcpy(__A, other.__A);
    }
    MoleculePart& operator=(const MoleculePart& other) {
        atoms_indices = other.atoms_indices;
        axis_indices  = other.axis_indices;
        gsl_vector_float_memcpy(__rot_axis, other.__rot_axis);
        gsl_vector_float_memcpy(__w, other.__w);
        gsl_matrix_float_memcpy(__A, other.__A);
        return *this;
    }

    void rotate_axis(const float angle, Atom * atom, gsl_vector_float * axis) {
        gsl_vector_float_set_zero(__w);
        set_quaternion_matrix(axis, angle, __A);
        for (auto j : atoms_indices) {
            rotate_vec((atom+j)->pos, __A, (atom+axis_indices.at(0))->pos, __w);
        };
        return;
    }
};

class AroundAxis : public MoleculePart {
    using MoleculePart::MoleculePart;
    public:
    void rotate(const float angle, Atom * atom) {
        gsl_vector_float_memcpy(__rot_axis, (atom+axis_indices.at(1))->pos);
        gsl_vector_float_sub(__rot_axis, (atom+axis_indices.at(0))->pos);
        gsl_vector_float_scale(__rot_axis, 1.0/gsl_blas_snrm2(__rot_axis));
        rotate_axis(angle, atom, __rot_axis);
    }
};

class FreeEnd : public MoleculePart {
    using MoleculePart::MoleculePart;
    public:
    void rotate(const float angle, Atom * atom, gsl_vector_float * axis) {
        rotate_axis(angle, atom, axis);
    }
};

class Molecule : public AtomHolder {
    public:
    using AtomHolder::AtomHolder;
    vector<AroundAxis> parts_around_axis;
    vector<FreeEnd> parts_free_end;
    // beware the constructor & copy assignment
    //gsl_vector_float * en_within;
};

class Surface: public AtomHolder {
    public:
    using AtomHolder::AtomHolder;
    // beware the constructor & copy assignment
    gsl_vector_float * v1;
    gsl_vector_float * v2;
};






void split_line(string line, string * splitted_line) {
    int pos=0, prev_state=0, curr_state=0;
    int i=0, j=0;
    for (auto cc: line) {
        curr_state = ((cc == ' ' || cc == '\t') ? 0 : 1);
        if (curr_state-prev_state == 1) {
            pos = i;
        } else if (curr_state-prev_state == -1) {
            *(splitted_line+j) = line.substr(pos, i-pos);
            ++j;
        }
        prev_state = curr_state;
        ++i;
    }
    if (prev_state == 1) {
        *(splitted_line+j) = line.substr(pos, i-pos);
    }
    return;
};
template <typename T> T * create_object_from_xyz(ifstream & file, string file_name) {
    string n_atoms;
    file.open(file_name);
    getline(file, n_atoms); // number of atoms
    T * molecule = new T(stoi(n_atoms));
    return molecule;
}
void load_xyz_to_object(ifstream & file, AtomHolder * holder) {
    string * splitted_line = new string[7];
    string line;
    int j = 0;
    while (getline(file, line)) {
        *(splitted_line+4) = "0.0";
        *(splitted_line+5) = "N";
        *(splitted_line+6) = "-1";
        split_line(line, splitted_line);
        string element = *splitted_line;
        float x = stof(*(splitted_line+1));
        float y = stof(*(splitted_line+2));
        float z = stof(*(splitted_line+3));
        float charge = stof(*(splitted_line+4));
        string hb_type = *(splitted_line+5);
        int hb_acceptor = stoi(*(splitted_line+6));
        (holder->atom+j)->fill(element, x,y,z, charge, hb_type, hb_acceptor);
        ++j;
    }
    file.close();
    delete[] splitted_line;
}
void load_atoms_from_xyz(ifstream & file, AtomHolder * holder, int n_atoms) {
    string * splitted_line = new string[7];
    string line;
    for (int j=0; j<n_atoms; ++j) {
        getline(file, line);
        *(splitted_line+4) = "0.0";
        *(splitted_line+5) = "N";
        *(splitted_line+6) = "-1";
        split_line(line, splitted_line);
        string element = *splitted_line;
        float x = stof(*(splitted_line+1));
        float y = stof(*(splitted_line+2));
        float z = stof(*(splitted_line+3));
        float charge = stof(*(splitted_line+4));
        string hb_type = *(splitted_line+5);
        int hb_acceptor = stoi(*(splitted_line+6));
        (holder->atom+j)->fill(element, x,y,z, charge, hb_type, hb_acceptor);
    }
    delete[] splitted_line;
}
void extract_basis(gsl_vector_float * v1, gsl_vector_float * v2, string title) {
    //std::cout << title.length() << endl;
    //std::cout << title << endl;
    char space = ' ';
    // 1 -- cursor inside number, 0 -- cursor outside number
    int bool_in_out_curr, bool_in_out_prev = 0;
    int num_begin_index, num_end_index;
    string num_str;
    float num;
    int n_nums = 0;
    for (int i=0; i<title.length(); ++i) {
        bool_in_out_curr = title[i] == space ? 0 : 1;
        //std::cout << bool_in_out_prev << ", " << bool_in_out_curr << ": ";
        if (bool_in_out_prev == 0 && bool_in_out_curr == 1) {
            // from space to digit
            num_begin_index = i;
            //std::cout << i;
        } else if (bool_in_out_prev == 1 && bool_in_out_curr == 0) {
            // from digit to space
            num_end_index = i;
            //std::cout << i;
            num_str = title.substr(num_begin_index, num_end_index-num_begin_index);
            num = std::stof(num_str);
            n_nums < 3 ? gsl_vector_float_set(v1, n_nums, num) : gsl_vector_float_set(v2, n_nums-3, num);
            //std::cout << " " << num;
            ++n_nums;
        }
        //std::cout << endl;
        bool_in_out_prev = bool_in_out_curr;
    }
    if (bool_in_out_prev == 1) {
        num_end_index = title.length();
        //std::cout << num_end_index;
        num_str = title.substr(num_begin_index, num_end_index-num_begin_index);
        num = std::stof(num_str);
        n_nums < 3 ? gsl_vector_float_set(v1, n_nums, num) : gsl_vector_float_set(v2, n_nums-3, num);
        //std::cout << " " << num << endl;
        ++n_nums;
    }
    cout << "v1 basis vector: ";
    print_gsl_vector(v1);
    cout << "v2 basis vector: ";
    print_gsl_vector(v2);
}
Surface * load_xyz_cell(string file_name, int read_basis) {
    ifstream file;
    Surface * surf = create_object_from_xyz<Surface>(file, file_name);
    string title;
    getline(file, title); // title
    surf->v1 = gsl_vector_float_calloc(3);
    surf->v2 = gsl_vector_float_calloc(3);
    if (read_basis == true) {
        // extract basis vectors from xyz file
        extract_basis(surf->v1, surf->v2, title);
    }
    load_xyz_to_object(file, surf);

    return surf;
}
Molecule * load_xyz_mol(string file_name) {
    ifstream file;
    Molecule * mol = create_object_from_xyz<Molecule>(file, file_name);
    string title;
    getline(file, title); // title
    load_xyz_to_object(file, mol);

    return mol;
}
list<int> extract_n_atoms(string title) {
    list<int> nums;
    //std::cout << title.length() << endl;
    //std::cout << title << endl;
    char space = ' ';
    // 1 -- cursor inside number, 0 -- cursor outside number
    int bool_in_out_curr, bool_in_out_prev = 0;
    int num_begin_index, num_end_index;
    string num_str;
    float num;
    int n_nums = 0;
    for (int i=0; i<title.length(); ++i) {
        bool_in_out_curr = title[i] == space ? 0 : 1;
        //std::cout << bool_in_out_prev << ", " << bool_in_out_curr << ": ";
        if (bool_in_out_prev == 0 && bool_in_out_curr == 1) {
            // from space to digit
            num_begin_index = i;
            //std::cout << i;
        } else if (bool_in_out_prev == 1 && bool_in_out_curr == 0) {
            // from digit to space
            num_end_index = i;
            //std::cout << i;
            num_str = title.substr(num_begin_index, num_end_index-num_begin_index);
            num = std::stoi(num_str);
            nums.push_back(num);
            //std::cout << " " << num;
            ++n_nums;
        }
        //std::cout << endl;
        bool_in_out_prev = bool_in_out_curr;
    }
    if (bool_in_out_prev == 1) {
        num_end_index = title.length();
        //std::cout << num_end_index;
        num_str = title.substr(num_begin_index, num_end_index-num_begin_index);
        num = std::stoi(num_str);
        nums.push_back(num);
        //std::cout << " " << num << endl;
        ++n_nums;
    }
    return nums;
}


void erase_file_contents(string path) {
    ofstream file;
    file.open(path);
    file.close();
    return;
    }
void open_file(ofstream & file, string type, string file_name) {
    if (type == "a") {
        file.open(file_name, std::ios_base::app);
    } else {
        file.open(file_name);
    }
    return;
}
void write_atoms_by_line(AtomHolder * molecule, ofstream & file) {
    Atom * atom;
    gsl_vector_float * pos;
    for (int i=0; i<molecule->n_atoms; ++i) {
        atom = molecule->atom+i;
        pos = atom->pos;
        file << atom->element << "\t";
        file << gsl_vector_float_get(pos, 0) << "\t";
        file << gsl_vector_float_get(pos, 1) << "\t";
        file << gsl_vector_float_get(pos, 2) << "\t";
        if (atom->hb_type == "N" && atom->hb_acceptor == -1) {
            file << atom->charge << "\n";
        } else {
            file << atom->charge << "\t";
            file << atom->hb_type << "\t";
            file << atom->hb_acceptor << "\n";
        }
    }
    return;
}
void write_xyz_single(AtomHolder * molecule, string file_name, string type="w", string description="--") {
    ofstream file;
    open_file(file, type, file_name);
    file << molecule->n_atoms << "\n";
    file << description << "\n";
    write_atoms_by_line(molecule, file);
    file.close();
    return;
}
void write_xyz_single_surface(Molecule * molecule, Surface * surface, string file_name, string type="w", string description="--") {
    ofstream file;
    open_file(file, type, file_name);
    file << molecule->n_atoms + surface->n_atoms << "\n";
    file << description << "\n";
    write_atoms_by_line(surface, file);
    write_atoms_by_line(molecule, file);
    file.close();
    return;
}
template <typename T> void write_xyz_list(list<T*> molecule_list, string file_name, string type="w", string description="--") {
    int n_atoms_tot=0;
    typename list<T*>::iterator iter = molecule_list.begin();
    for (iter; iter!=molecule_list.end(); ++iter) {
        n_atoms_tot += (*iter)->n_atoms;
    }
    ofstream file;
    open_file(file, type, file_name);
    file << n_atoms_tot << "\n";
    file << description << "\n";
    int j = 0;
    iter = molecule_list.begin();
    for (iter; iter!=molecule_list.end(); ++iter, ++j) {
        write_atoms_by_line(*iter, file);
    }
    file.close();
    return;
}
void write_xyz_list_surface(list<Molecule*> molecule_list, Surface * surface, string file_name, string type="w", string description="--") {
    int n_atoms_tot=surface->n_atoms;
    list<Molecule*>::iterator iter = molecule_list.begin();
    for (iter; iter!=molecule_list.end(); ++iter) {
        n_atoms_tot += (*iter)->n_atoms;
    }
    ofstream file;
    open_file(file, type, file_name);
    file << n_atoms_tot << "\n";
    file << description << "\n";
    write_atoms_by_line(surface, file);
    int j = 0;
    iter = molecule_list.begin();
    for (iter; iter!=molecule_list.end(); ++iter, ++j) {
        write_atoms_by_line(*iter, file);
    }
    file.close();
    return;
}


void write_atoms_by_line_offset(AtomHolder * molecule, ofstream & file, gsl_vector_float * offset) {
    Atom * atom;
    gsl_vector_float * v = gsl_vector_float_alloc(3);
    for (int i=0; i<molecule->n_atoms; ++i) {
        atom = molecule->atom+i;
        gsl_vector_float_memcpy(v, offset);
        gsl_vector_float_add(v, atom->pos);
        file << atom->element << "\t";
        file << gsl_vector_float_get(v, 0) << "\t";
        file << gsl_vector_float_get(v, 1) << "\t";
        file << gsl_vector_float_get(v, 2) << "\t";
        if (atom->hb_type == "N" && atom->hb_acceptor == -1) {
            file << atom->charge << "\n";
        } else {
            file << atom->charge << "\t";
            file << atom->hb_type << "\t";
            file << atom->hb_acceptor << "\n";
        }
    }
    gsl_vector_float_free(v);
    return;
}

void write_xyz_single_periodic(gsl_vector_float * v1, gsl_vector_float * v2, AtomHolder * surface, string file_name, string type="w", string description="--") {
    ofstream file;
    open_file(file, type, file_name);
    file << 9*(surface->n_atoms) << "\n";
    file << description << "\n";
    vector<int> mults = {1,0, 1,1, 0,1, -1,1, -1,0, -1,-1, 0,-1, 1,-1};
    gsl_vector_float * w = gsl_vector_float_alloc(3);
    write_atoms_by_line(surface, file);
    for (int i=0; i<(mults.size()/2); ++i) {
        gsl_vector_float_set_zero(w);
        gsl_blas_saxpy(mults.at(2*i), v1, w);
        gsl_blas_saxpy(mults.at(2*i+1), v2, w);
        write_atoms_by_line_offset(surface, file, w);
    }
    file.close();
    gsl_vector_float_free(w);
    return;
}
void write_xyz_single_surface_periodic(gsl_vector_float * v1, gsl_vector_float * v2, Molecule * molecule, Surface * surface, string file_name, string type="w", string description="--") {
    ofstream file;
    open_file(file, type, file_name);
    file << 9*(molecule->n_atoms + surface->n_atoms) << "\n";
    file << description << "\n";
    vector<int> mults = {1,0, 1,1, 0,1, -1,1, -1,0, -1,-1, 0,-1, 1,-1};
    gsl_vector_float * w = gsl_vector_float_alloc(3);
    write_atoms_by_line(molecule, file);
    write_atoms_by_line(surface, file);
    for (int i=0; i<(mults.size()/2); ++i) {
        gsl_vector_float_set_zero(w);
        gsl_blas_saxpy(mults.at(2*i), v1, w);
        gsl_blas_saxpy(mults.at(2*i+1), v2, w);
        write_atoms_by_line_offset(molecule, file, w);
        write_atoms_by_line_offset(surface, file, w);
    }
    file.close();
    gsl_vector_float_free(w);
    return;
}
void write_xyz_list_surface_periodic(gsl_vector_float * v1, gsl_vector_float * v2, list<Molecule*> molecule_list, Surface * surface, string file_name, string type="w", string description="--") {
    int n_atoms_tot=surface->n_atoms;
    list<Molecule*>::iterator iter = molecule_list.begin();
    for (iter=molecule_list.begin(); iter!=molecule_list.end(); ++iter) {
        n_atoms_tot += (*iter)->n_atoms;
    }
    ofstream file;
    open_file(file, type, file_name);
    file << 9*n_atoms_tot << "\n";
    file << description << "\n";
    vector<int> mults = {0,0, 1,0, 1,1, 0,1, -1,1, -1,0, -1,-1, 0,-1, 1,-1};
    gsl_vector_float * w = gsl_vector_float_alloc(3);
    for (int i=0; i<(mults.size()/2); ++i) {
        gsl_vector_float_set_zero(w);
        gsl_blas_saxpy(mults.at(2*i), v1, w);
        gsl_blas_saxpy(mults.at(2*i+1), v2, w);
        write_atoms_by_line_offset(surface, file, w);
        for (iter=molecule_list.begin(); iter!=molecule_list.end(); ++iter) {
            write_atoms_by_line_offset(*iter, file, w);
        }
    }
    file.close();
    gsl_vector_float_free(w);
    return;
}

void reverse_mol(Molecule * mol) {

    for (int i=0; i<mol->n_atoms; ++i) {
        //print_gsl_vector((mol->atom+i)->pos);
        gsl_vector_float_set((mol->atom+i)->pos,1,-1*gsl_vector_float_get((mol->atom+i)->pos,1));
        //print_gsl_vector((mol->atom+i)->pos);
    }

    return;
}