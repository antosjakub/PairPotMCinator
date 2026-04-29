#include <fstream>
#include <iostream>
#include <string>
#include <filesystem>
#include "finder.h"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

void write_blank_line(ofstream file) {
    file << std::endl;
}

template <typename T>
map<string, T> DOF_base {
    {"translate_2D", 0},
    {"translate_e3", 0},
    {"rotate_3D", 0},
    {"rotate_e3", 0},
    {"deform", 0}
};



int main(int argc, char * argv[]) {
    int DEBUG = 1;
    consts::init();

    // 0) READ INPUT 
    ifstream data_file(argv[1]);
    json data_json = json::parse(data_file);

    // Check that the .json file contains all the necessary high-level sections.
    auto raise_missing_section = [](string section_name) {
        std::cout << "ERROR: You forgot to include the " << section_name << " section in your .json file." << std::endl;
    };
    if (data_json["unit_cell"].is_null()) {
        raise_missing_section("'unit_cell'");
    };
    if (data_json["molecules"].is_null() && data_json["molecule_list"].is_null()) {
        raise_missing_section("'molecules' or 'molecule_list'");
    };
    if (data_json["results"].is_null()) {
        raise_missing_section("'results'");
    };

    // 1) CONSTRUCT THE SURFACE FROM THE UNIT CELL
    std::cout << "====== Reading cell data..." << std::endl;
    auto cell_path = data_json["unit_cell"]["path"].get<std::string>();
    Surface * cell_1x1 = load_xyz_cell(cell_path, true);
    cell_1x1->center();
    Surface * cell;
    if (!data_json["unit_cell"]["multiplicity"].is_null()) {
        auto cell_multiplicity = data_json["unit_cell"]["multiplicity"].get<std::vector<int>>();
        cell = return_cell_multiple(cell_1x1, cell_multiplicity.at(0), cell_multiplicity.at(1));
        cell->center();
    } else {
        cell = cell_1x1;
    }
    std::cout << "====== Cell data read sucessfully." << std::endl;


    // 2) INIT NECESSARY VARIABLES
    list<Molecule*> mol_list;
    json mol_dict;


    // 3) CONSTRUCT MOLECULES
    // True if "molecules" key present
    bool bool_molecules = !data_json["molecules"].is_null();
    // True if "molecule_list" key present
    bool bool_molecules_list = !data_json["molecule_list"].is_null();
    if (bool_molecules and not bool_molecules_list) {
        std::cout << "====== Prepating list of molecules..." << std::endl;
        // Goal: do not compute descent mol multiple times 
        // a) load mols - no copies, declare deform, etc.
        // b) descent them one by one
        // c) create mol_list <- contains copies - full
        std::cout << "==== Reading molecules from 'molecules' key..." << std::endl;
        Molecule * mol;
        list<Molecule*> mol_list_no_copies;
        mol_dict = data_json["molecules"];
        for (auto& mol_info : data_json["molecules"]) {
            std::cout << "molecule_name: " << mol_info["name"] << std::endl;
            mol = load_xyz_mol(mol_info["path"]);
            mol->center();
            if (!mol_info["deformation"].is_null()) {
                std::cout << "use deform:" << std::endl;
                for (auto & deform_part : mol_info["deformation"]) {
                    std::cout << deform_part["type"] << std::endl;
                    auto type = deform_part["type"].get<std::string>();
                    auto axis_indices = deform_part["axis_indices"].get<std::vector<int>>();
                    auto atom_indices = deform_part["atom_indices"].get<std::vector<int>>();
                    substract_one_from_vector(axis_indices);
                    substract_one_from_vector(atom_indices);
                    if (type == "around_axis") {
                        mol->parts_around_axis.push_back(AroundAxis(axis_indices, atom_indices));
                    } else if (type == "free_end") {
                        mol->parts_free_end.push_back(FreeEnd(axis_indices, atom_indices));
                    } else {
                        std::cout << "WARNING: Incorrect deform type passed" << std::endl;
                    }
                }
            }
            mol_list_no_copies.push_back(mol);
        }
        std::cout << "==== Molecules read sucessfully." << std::endl;

        std::cout << std::endl;
        Options::init_mol_list_surf(mol_list_no_copies, cell);
        PairPotentials::init();

        std::cout << "==== Descending molecules..." << std::endl;
        list<Molecule*>::iterator iter;
        for (iter=mol_list_no_copies.begin(); iter!=mol_list_no_copies.end(); ++iter) {
            descent_mol(*iter, cell);
        }
        std::cout << "==== Descended molecules sucessfully." << std::endl;

        std::cout << "==== Creating full list of molecules..." << std::endl;
        iter=mol_list_no_copies.begin();
        for (auto& mol_info : data_json["molecules"]) {
            append_mols_to_mol_list(mol_list, *iter, mol_info["n_copies"]);
            ++iter;
        }
        std::cout << "==== Created full list of molecules." << std::endl;

        std::cout << "====== List of molecules prepared." << std::endl;


    } else if (bool_molecules_list and not bool_molecules) {

        std::cout << "====== Prepating list of molecules..." << std::endl;
        // a) read mol list from xyz
        // b) iterate over the list of atom counts (ex. "23 34 21")
        // c) for each count: create new mol, read lines & populate atoms -> mol_list

        std::cout << "==== Reading molecule atom counts from .xyz file..." << std::endl;
        ifstream file;
        file.open(data_json["molecule_list"]["path"]);
        string n_atoms_total, n_atoms_string;
        getline(file, n_atoms_total);  // 1st line
        getline(file, n_atoms_string); // 2nd line
        list<int> n_atoms_list = extract_n_atoms(n_atoms_string);
        std::cout << "==== Molecule atom counts read." << std::endl;

        std::cout << "==== Loading molecule objects from .xyz file..." << std::endl;
        list<int>::iterator iter_atom_count = n_atoms_list.begin();
        for (iter_atom_count; iter_atom_count!=n_atoms_list.end(); ++iter_atom_count) {
            Molecule * mol = new Molecule(*iter_atom_count);
            load_atoms_from_xyz(file, mol, *iter_atom_count);
            mol->set_center_of_mass(mol->pos);
            std::cout << "pos of mol: ";
            print_gsl_vector(mol->pos);
            mol_list.push_back(mol);
        }
        file.close();
        std::cout << "==== Molecule objects from .xyz file loaded." << std::endl;

        std::cout << "==== Populating molecule objects with data specified in the json file..." << std::endl;
        std::cout << std::endl;
        std::cout << "Reading data from json..." << std::endl;
        list<Molecule*>::iterator iter_mol = mol_list.begin();
        list<Molecule*>::iterator iter_mol_base = mol_list.begin();
        mol_dict = data_json["molecule_list"]["molecules"];
        for (auto& mol_info : data_json["molecule_list"]["molecules"]) {
            std::cout << "molecule_name: " << mol_info["name"] << std::endl;
            if (!mol_info["deformation"].is_null()) {
                std::cout << "- deformation declared:" << std::endl;
                for (auto & deform_part : mol_info["deformation"]) {
                    iter_mol = iter_mol_base;
                    std::cout << deform_part["type"] << std::endl;
                    auto type = deform_part["type"].get<std::string>();
                    auto axis_indices = deform_part["axis_indices"].get<std::vector<int>>();
                    auto atom_indices = deform_part["atom_indices"].get<std::vector<int>>();
                    substract_one_from_vector(axis_indices);
                    substract_one_from_vector(atom_indices);
                    if (type == "around_axis") {
                        for (int i=0; i<mol_info["n_copies"]; ++i) {
                            (*iter_mol)->parts_around_axis.push_back(AroundAxis(axis_indices, atom_indices));
                            ++iter_mol;
                        };
                    } else if (type == "free_end") {
                        for (int i=0; i<mol_info["n_copies"]; ++i) {
                            (*iter_mol)->parts_free_end.push_back(FreeEnd(axis_indices, atom_indices));
                            ++iter_mol;
                        }
                    } else {
                        std::cout << "WARNING: Incorrect deform type passed" << std::endl;
                    }
                }
            } else {
                iter_mol = iter_mol_base;
                for (int i=0; i<mol_info["n_copies"]; ++i) {
                    ++iter_mol;
                }
            }
            iter_mol_base = iter_mol;
        }
        std::cout << "==== Molecule objects populated with json data." << std::endl;

        std::cout << std::endl;
        Options::init_mol_list_surf(mol_list, cell);
        PairPotentials::init();
    }
    std::cout << std::endl;
    

    auto n_cycles = data_json["n_cycles"].get<int>();
    // scatter ?
    int use_scatter = false;
    Stage stage_scatter("scatter", 1, 42); // create object outside of this if because otherwise the if later on is going to complain
    if (!data_json["scatter"].is_null()) {
        use_scatter = true;
        stage_scatter.n_trials = data_json["scatter"]["n_trials"];
        map<string, float> universal_params = DOF_base<float>;
        if (data_json["scatter"]["translate_2D"] == true) {
            float v1rand = 0.5*gsl_blas_snrm2(cell->v1);
            float v2rand = 0.5*gsl_blas_snrm2(cell->v2);
            universal_params["translate_2D"] = 0.5*(v1rand + v2rand);
        }
        if (data_json["scatter"]["translate_e3"] == true) {
            universal_params["translate_e3"] = 1.0;
        }
        if (data_json["scatter"]["rotate_e3"] == true) {
            universal_params["rotate_e3"] = 1.0;
        }
        if (data_json["scatter"]["rotate_3D"] == true) {
            universal_params["rotate_3D"] = 1.0;
        }
        if (data_json["scatter"]["deform"] == true) {
            universal_params["deform"] = 1.0;
        }
        // iterate over all defined molecule groups and create DOF dict for bind and params
        for (auto& mol_info : mol_dict) {
            map<string, float> new_group_params = universal_params;
            stage_scatter.params.insert(make_pair(mol_info["name"], new_group_params));
        }
        // lets print out the map
        std::cout << std::endl;
        std::cout << "SCATTER" << std::endl;
        for (auto group : stage_scatter.params) {
            std::cout << "Molecule: " << group.first << std::endl;
            for (auto dof : stage_scatter.params[group.first]) {
                std::cout << "  " << dof.first << ": " << dof.second << std::endl;
            }
        }
        //stage_scatter.print_out();
    } else {
        //delete[] &stage_scatter;
    }
    std::cout << "use_scatter: " << use_scatter << std::endl;

    // DECALARE STAGES
    vector<Stage> stages;
    for (auto& stage_js : data_json["stages"]) {
        //std::cout << std::endl;
        //std::cout << stage_js["params"] << std::endl;

        Stage stage(stage_js["name"], stage_js["n_loops"], stage_js["n_trials"]);

        map<string, float> universal_params = DOF_base<float>;

        if (!stage_js["params"].is_null()) {
            if (!stage_js["params"]["translate_e3"].is_null()) {
                universal_params["translate_e3"] = stage_js["params"]["translate_e3"];
            }
            if (!stage_js["params"]["rotate_e3"].is_null()) {
                universal_params["rotate_e3"] = stage_js["params"]["rotate_e3"];
            }
            if (!stage_js["params"]["rotate_3D"].is_null()) {
                universal_params["rotate_3D"] = stage_js["params"]["rotate_3D"];
            }
            if (!stage_js["params"]["deform"].is_null()) {
                universal_params["deform"] = stage_js["params"]["deform"];
            }
            if (!stage_js["params"]["translate_2D"].is_null()) {
                universal_params["translate_2D"] = stage_js["params"]["translate_2D"];
            }
        }

        // iterate over all defined molecule groups and create DOF dict for bind and params
        for (auto& mol_info : mol_dict) {
            string mol_name = mol_info["name"];
            map<string, float> new_group_params = universal_params;
            stage.params.insert(make_pair(mol_name, new_group_params));
            map<string, int> new_group_bind = DOF_base<int>;
            stage.bind.insert(make_pair(mol_name, new_group_bind));
            // is the molecule group name in the "params" section? 
            for (auto& [param_group, params] : stage_js["params"].items()) {
                if (param_group == mol_name) {
                    if (!stage_js["params"][param_group]["translate_e3"].is_null()) {
                        stage.params[param_group]["translate_e3"] = stage_js["params"][param_group]["translate_e3"];
                    }
                    if (!stage_js["params"][param_group]["rotate_e3"].is_null()) {
                        stage.params[param_group]["rotate_e3"] = stage_js["params"][param_group]["rotate_e3"];
                    }
                    if (!stage_js["params"][param_group]["rotate_3D"].is_null()) {
                        stage.params[param_group]["rotate_3D"] = stage_js["params"][param_group]["rotate_3D"];
                    }
                    if (!stage_js["params"][param_group]["deform"].is_null()) {
                        stage.params[param_group]["deform"] = stage_js["params"][param_group]["deform"];
                    }
                    if (!stage_js["params"][param_group]["translate_2D"].is_null()) {
                        stage.params[param_group]["translate_2D"] = stage_js["params"][param_group]["translate_2D"];
                    }
                }
            }
        }
        if (!stage_js["bind"].is_null()) {
            // iterate over molecule groups mentioned in "bind"
            for (auto& [group_name, DOF_list] : stage_js["bind"].items()) {
                auto DOF_vect = DOF_list.get<vector<string>>();
                for (auto dof : DOF_vect) {
                    stage.bind[group_name][dof] = 1;
                }
            }
        }
        if (DEBUG == 1) {
            std::cout << std::endl;
            std::cout << "NEW STAGE: " << stage.name << std::endl;
            std::cout << "======== START: Now printing params." << std::endl;
            for (auto group : stage.params) {
                std::cout << "Molecule: " << group.first << std::endl;
                for (auto dof : stage.params[group.first]) {
                    std::cout << "  " << dof.first << ": " << dof.second << std::endl;
                }
            }
            std::cout << "======== END: Now printing params." << std::endl;
            std::cout << "======== START: Now printing bind." << std::endl;
            for (auto group : stage.bind) {
                std::cout << "Molecule: " << group.first << std::endl;
                for (auto dof : stage.bind[group.first]) {
                    std::cout << "  " << dof.first << ": " << dof.second << std::endl;
                }
            }
            std::cout << "======== END: Now printing bind." << std::endl;
        }
        std::cout << std::endl;
        //stage.print_out();
        stages.push_back(stage);
    }



    // PREPARE RESULTS FOLDER AND FILES 
    string result_folder = data_json["results"]["folder"];
    std::filesystem::create_directory(result_folder);
    map<string, map<string, LogObject*>> log;
    for (const auto& [log_type, log_list] : log_source) {
        string result_folder_type = result_folder + log_type+"_files/";
        std::filesystem::create_directory(result_folder_type);
        for (auto const & log_var_name : log_list) {
            string log_var_path = result_folder_type + log_var_name + "." + log_type;
            bool log_var_use = data_json["results"]["."+log_type][log_var_name];
            if (log_var_use) erase_file_contents(log_var_path);
            log[log_type][log_var_name] = new LogObject(log_var_use, log_var_path);
        }
    };


    // {"DTDPP": 2, "DPDPP": 4} <-- number of molecules in the group
    map<string, int> group_counts;
    for (auto& mol_info : mol_dict) {
        group_counts.insert(make_pair(mol_info["name"], mol_info["n_copies"]));
    };

    // FINDER
    auto rand_seed = data_json["random_seed"].get<int>();
    Finder finder(cell, mol_list, group_counts, n_cycles, rand_seed, log);

    string string_cycle;
    gsl_vector_float * en_cycles = gsl_vector_float_calloc(n_cycles);
    if (log["txt"]["cycle_min"]->use)  open_file(finder.log["txt"]["cycle_min"]->file, "w", log["txt"]["cycle_min"]->path);
    if (log["txt"]["stages_min"]->use) open_file(finder.log["txt"]["stages_min"]->file, "w", log["txt"]["stages_min"]->path);
    if (log["txt"]["all_min"]->use)    open_file(finder.log["txt"]["all_min"]->file, "w", log["txt"]["all_min"]->path);

    if (log["txt"]["cycle_min"]->use)  finder.log["txt"]["cycle_min"]->file << std::endl;
    for (int i=0; i<n_cycles; ++i) {
        // Log.
        string_cycle = "cycle: " + to_string(i+1);
        std::cout << std::endl;
        if (log["txt"]["all_min"]->use) {
            finder.log["txt"]["all_min"]->file << std::endl;
        }
        finder.prepare_for_new_cycle();

        if (log["txt"]["stages_min"]->use) finder.log["txt"]["stages_min"]->file << std::endl;

        if (use_scatter) {
            if (bool_molecules) {
                finder.run_scatter(stage_scatter, string_cycle + ", scatter");
                if (log["txt"]["stages_min"]->use) finder.log["txt"]["stages_min"]->file << string_cycle + ", scatter" << " | en: " << finder.en_min << std::endl;
                if (log["xyz"]["stages_min_cell"]->use) write_xyz_list_surface(finder.mol_list_curr, cell, log["xyz"]["stages_min_cell"]->path, "a", string_cycle + ", stage: scatter");
                if (log["xyz"]["stages_min_periodic"]->use) write_xyz_list_surface_periodic(cell->v1, cell->v2, finder.mol_list_curr, cell, log["xyz"]["stages_min_periodic"]->path, "a", string_cycle + ", stage: scatter");
            } else {
                std::cout << "WARNING: you are using scatter even though a list of molecules was given as an input. The molecules will be initialized randomly.";
            }
        } else {
            if (bool_molecules_list) {
                finder.en_min = finder.calc_en_mol_list_periodic(mol_list);
            } else {
                finder.en_min = std::numeric_limits<int>::max();
            }
            std::cout << string_cycle + ", init" << " | en: " << finder.en_min << std::endl;
            if (log["txt"]["all_min"]->use) finder.log["txt"]["all_min"]->file << string_cycle + ", init" << " | en: " << finder.en_min << std::endl;
            if (log["xyz"]["all_min_cell"]->use) write_xyz_list_surface(finder.mol_list_curr, cell, log["xyz"]["all_min_cell"]->path, "a", string_cycle+", init");
            if (log["xyz"]["all_min_periodic"]->use) write_xyz_list_surface_periodic(cell->v1, cell->v2, finder.mol_list_curr, cell, log["xyz"]["all_min_periodic"]->path, "a", string_cycle+", init");
            if (log["txt"]["stages_min"]->use) finder.log["txt"]["stages_min"]->file << string_cycle + ", init" << " | en: " << finder.en_min << std::endl;
            if (log["xyz"]["stages_min_cell"]->use) write_xyz_list_surface(finder.mol_list_curr, cell, log["xyz"]["stages_min_cell"]->path, "a", string_cycle + ", stage: init");
            if (log["xyz"]["stages_min_periodic"]->use) write_xyz_list_surface_periodic(cell->v1, cell->v2, finder.mol_list_curr, cell, log["xyz"]["stages_min_periodic"]->path, "a", string_cycle + ", stage: init");
        }

        // STAGES
        for(int j=0; j<stages.size(); ++j) {
            std::cout << string_cycle << "  -----------------------" << std::endl;
            if (log["txt"]["all_min"]->use) finder.log["txt"]["all_min"]->file << string_cycle << "  -----------------------" << std::endl;

            finder.run_stage(stages.at(j), string_cycle + ", stage: " + to_string(j+1));

            if (log["txt"]["stages_min"]->use) finder.log["txt"]["stages_min"]->file << string_cycle + ", stage: " + to_string(j+1) << " | en: " << finder.en_min << std::endl;
            if (log["xyz"]["stages_min_cell"]->use) write_xyz_list_surface(finder.mol_list_curr, cell, log["xyz"]["stages_min_cell"]->path, "a", "cycle: " + to_string(i+1) + ", stage: " + to_string(j+1));
            if (log["xyz"]["stages_min_periodic"]->use) write_xyz_list_surface_periodic(cell->v1, cell->v2, finder.mol_list_curr, cell, log["xyz"]["stages_min_periodic"]->path, "a", "cycle: " + to_string(i+1) + ", stage: " + to_string(j+1));
        }


        // RESULTS FOR CYCLE
        std::cout << string_cycle << "  -----------------------" << std::endl;
        std::cout << string_cycle + ", cycle_min: " << finder.en_min << std::endl;
        gsl_vector_float_set(en_cycles, i, finder.en_min);
        if (log["txt"]["all_min"]->use) finder.log["txt"]["all_min"]->file << string_cycle << "  -----------------------" << std::endl;
        if (log["txt"]["all_min"]->use) finder.log["txt"]["all_min"]->file << string_cycle + ", cycle_min: " << finder.en_min << std::endl;
        if (log["txt"]["cycle_min"]->use) finder.log["txt"]["cycle_min"]->file << string_cycle + " | en: " << finder.en_min << std::endl;
        if (log["xyz"]["cycle_min_cell"]->use) write_xyz_list_surface(finder.mol_list_curr, cell, log["xyz"]["cycle_min_cell"]->path, "a", "cycle: " + to_string(i+1));
        if (log["xyz"]["cycle_min_periodic"]->use) write_xyz_list_surface_periodic(finder.v1, finder.v2, finder.mol_list_curr, cell, log["xyz"]["cycle_min_periodic"]->path, "a", "cycle: " + to_string(i+1));
    }


    return 0;

}