


PairPotMCinator: a tool for fast simulations of the organic-molecular ordering on the solid surfaces using pair potentials

# Supplementary


### TD
- systematic benchmarking - report both success and failure
- criteria for convergence
- energy diffs between competing configurations
- compare independent MC runs


show diff MC runs
- diff cycles? similar minima
- diff cycles? diff minima - hard energy landscape

konvergence záleží na množství stupňů volnosti a na nastavení cyklů, což je ale celkem zjevný. Zkusil bych na jednom/dvou příkladech ukázat konvergenci energie, jednou, kdy to ve většině cyklů dokonverguje do stejného minima, a jednou, kdy to je příliš složité a konverguje to vždy někam jinam. To, co generuješ v all_min.txt by mělo stačit...

### with FFT
for large number of runs - use FFT to extract only the periodic superstructures
Do 5 runs, each with 10 cycles
Do 10 runs, each with 100 cycles


### outline


criterium for convergence?
- enough configurations share the name structure

show real run - dpp 2x1 config, dpp 3x1 config
- energy diffs betweeen configuraitons
[text](results/DPP_on_graphene/config_DTDPP_1x3.json)
--- 

## trashbin

In all three cases, we consider 10 cycles. first n times independently scatter the molecules on the surface, select the configuration with minimal energy. We then proceed in two stages.


## Algorithm convergence

We illustrate the process of determining the convergence of the algorithm on the previously mentioned DPDPP molecules on graphene.

Consider the cells 1x1, 1x2, and 1x3.

### cell 1x2

In the case of the 1x2 unit cell, the algorithm arrives at the same molecular superstructure in the majority of cycles.

-- show the prevelant superstruct & energies for it ---
-- show all the energies of the configuration --
-- show some other ones

### cell 1x3

- multiple different structures - same energy
- 2 similar ones - minimum?
- run more, see which repeat

In this case, each MC step has more degrees of freedom that in the case of the 1x2 unit cell, hence we will need more runs to find similar structures


### summary

In short, We have convergence, if multiple superstructures have the same features

Manually go through. Or use FFT to detect similar ones, cluster them, ...

In case of more cycles - 100 we can you fft or randon transform to cluster similar structures
the repository provides both approaches in the post processing step


--------------------------------------------
--------------------------------------------
--------------------------------------------

### Algorithm convergence

We illustrate the process of determining the convergence of the algorithm on the previously mentioned DPDPP molecules on graphene with the unit cell of size 1x3.

We run the MC algorithm 5 times, each run with the following config file:

```
// config.json
"molecules": [{
    "name": "DTDPP",
    "path": "xyz_files/DTDPP.xyz",
    "n_copies": 3,
    "deformation": [{
            "name": "deform_1",
            "type": "free_end",
            "axis_indices": [5,7],
            "atom_indices": [7,8,9,10,11,12,13,14]
        }, {
            "name": "deform_2",
            "type": "free_end",
            "axis_indices": [19,21], 
            "atom_indices": [21,22,23,24,25,26,27,28]
        }
    ]
}],
"n_cycles": 10,
"scatter": {
    "n_trials": 64,
    "translate_2D": true,
    "translate_e3": false,
    "rotate_e3": true,
    "rotate_3D": false,
    "deform": false
},
"stages": [{
    "name": "stage_1",
    "n_loops": 200,
    "n_trials": 16,
    "params": {
        "translate_2D": 1.0,
        "rotate_e3": 0.1
    }
}, {
    "name": "stage_2",
    "n_loops": 400,
    "n_trials": 16,
    "params": {
        "translate_2D": 0.5,
        "translate_e3": 0.1,
        "rotate_3D": 0.05,
        "deform": 0.05
    }
}]
```

The table below shows the number of degrees of freedom for each operation on a single molecule of DTDPP:

operation, description, number of degrees of freedom
- translate_2D, translation in the plane parallel with the surface, 2
- translate_3e, translation in the direction perpendicular to the surface, 1
- rotate_3D, rotation about a random axis, 3
- deform, 2 rotations of part of the molecule about a random axis, 6 (2x3)

Together:
- stage 1: 3 mols, each with 2+1=3 degrees of freedom
- stage 2: 3 mols, each with 2+1+3+6=12 degrees of freedom


Combining the results from the five runs, we get 50 (5x10) configurations - corresponding to the final configurations of each cycle.

(We get a similar energy landscape, but different structures.)

-- show a plot
-- x-axis: index of the cycle
-- y-axis: optimal energy 

In case we simulate only a few cycles, we can go through the results manually with an aid of some molecular visualization software (such as jmol).
If we get similar structures across multiple cycles, we have converged into a stable minimum. If this is not case, we have to tweak the hyperparmeters of config file.

However, if we run the simulaiton for multiple cycles, going through the results manually is very timeconsuming. To aid in this endevour, we can use FFT, radon transform, or ML image classificaiton method to automatically cluster the results.

Our repository provides a few different tool in this direction.

Using, for example, the radon tranform, we get 5 different cluster correspoding to 5 main configuration patterns.

-- show a plot with examples of 3 configs --

To show an example where the algorithm struggle to converge into a some stable minima, we consider the following config file:

```
"n_cycles": 10,
"scatter": {
    "n_trials": 64,
    "translate_2D": true,
    "translate_e3": false,
    "rotate_e3": true,
    "rotate_3D": false,
    "deform": false
},
"stages": [{
    "name": "stage_1",
    "n_loops": 500,
    "n_trials": 16,
    "params": {
        "translate_2D": 0.5,
        "translate_e3": 0.1,
        "rotate_3D": 0.5,
        "deform": 0.5
    }
}]
```

Here, directly start with a stage containding defomation of the moleucles.
The simulation produces various different superstructes, without any significant between them.

-- show example of 4 configs --





passthrough A:

0: random
1: shoulder to shoulder
2: full col full row
3: fork

cycle 1:
2323133213

cycle 2:
2323133213

cycle 3:
2323032213

cycle 4:
3303132013

cycle 5:
2333132013




Now, onsider a case where, apart from a translation and rotation in a plane, we also allow the molecules to move and rotate in 3d space, and internally deform - all within a single stage:

\begin{verbatim}
"stages": [{
    "name": "planar",
    "n_loops": 200,
    "n_trials": 16,
    "params": {
        "translate_2D": 1.0,
        "rotate_e3": 0.1
    }
}, {
    "name": "planar_and_deform",
    "n_loops": 400,
    "n_trials": 16,
    "params": {
        "translate_2D": 0.5,
        "translate_e3": 0.1,
        "rotate_e3": 0.1,
        "rotate_3D": 0.1,
        "deform": 0.1
    }
}]
\end{verbatim}

We thus get
\begin{itemize}
    \item stage 1: 3 mols, each with $2+1+1+3+2\cdot3$ = 13 DoF, giving a total of 39 DoF
\end{itemize}

Again, we, get 10 optimal configurations. However this time, the molecules are all tangled up and even 4 more runs do not help us. See Fig. () for an example of a few such configurations.



-------- LABEL = 0 --------
1: run=3, cycle=5
2: run=4, cycle=3
3: run=4, cycle=8
4: run=5, cycle=8
---------------------------
-------- LABEL = 1 --------
1: run=1, cycle=5
2: run=1, cycle=9
3: run=2, cycle=5
4: run=2, cycle=9
5: run=3, cycle=9
6: run=4, cycle=5
7: run=4, cycle=9
8: run=5, cycle=5
9: run=5, cycle=9
---------------------------