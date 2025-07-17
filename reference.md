# Structure of the .json file

The .json file should contain the following 4 sections:


## Surface unit cell

By default:

```
"unit_cell": {
    "path": "xyz_files/cell.xyz"
},
```

If you want to periodically copy the given surface cell, for example, 2 times in the direction of the x-axis and 3 times in the direction of the y-axis:

```
"unit_cell": {
    "path": "xyz_files/cell.xyz",
    "multiplicity": [2, 3]
},
```


## Molecules

If you have a .xyz file for every molecule (use if you want to position the molecules randomly on the surface):

```
"molecules": [{
    "name": "DTDPP",
    "path": "xyz_files/DTDPP.xyz",
    "n_copies": 3
}, {
    "name": "DPDPP",
    "path": "xyz_files/DPDPP.xyz",
    "n_copies": 2
}],
```

If you have a the molecules stored in one .xyz (use if you want to initialize the molecules at specific positions):

```
"molecule_list": {
    "path": "xyz_files/mol_list.xyz",
    "molecules": [{
        "name": "DTDPP",
        "n_copies": 3,
    }, {
        "name": "DPDPP",
        "n_copies": 2,
    }]
},

```


## Simulation parameters

Parameters governing the Monte Carlo simulation:

```
"random_seed": 1,
"n_cycles": 2,
"scatter": {
    "n_trials": 32,
    "translate_2D": true,
    "translate_e3": false,
    "rotate_e3": true,
    "rotate_3D": false,
    "deform": false
},
"stages": [{
    "name": "planar",
    "n_loops": 100,
    "n_trials": 4,
    "params": {
        "translate_2D": 0.25,
        "rotate_e3": 0.25
    }
}, {
    "name": "settle_into_place",
    "n_loops": 50,
    "n_trials": 8,
    "params": {
        "translate_2D": 0.1,
        "rotate_e3": 0.1,
        "translate_e3": 0.1,
        "rotate_3D": 0.1
    }
}],
```


## Saving results

A section specifying at what level of detail to store the results.
The .xyz files store the position of each atom of the slab, the .txt files record the energy of the slab.

- 'cycle_min' = store the best individual from each cycle

- 'stage_min' = store the best individual from each stage

- 'all_min' = store every time a better individual is found


```
"results": {
    "folder": "results_simple/",
    ".xyz": {
        "cycle_min_periodic": true,
        "cycle_min_cell": true,
        "stages_min_periodic": true,
        "stages_min_cell": true,
        "all_min_periodic": true,
        "all_min_cell": true
    },
    ".txt": {
        "cycle_min": true,
        "stages_min": true,
        "all_min": true
    }
}
```
