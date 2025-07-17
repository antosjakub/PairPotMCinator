# PairPotMCinator

### Table of Contents

- [Dependencies](#dependencies)
- [Installing this project](#installing-this-project)
- [Example walkthrough](#example-walkthrough)
- [Editing the source code](#editing-the-source-code)


## Dependencies

### nlohmann
-- A C++ library providing I/O interface for .json files.
<details>
<summary>Installation steps</summary>
<br>

```
sudo apt-get install nlohmann-json3-dev
```
</details>

### Cmake, g++, make
-- A set of linux tools used for building source code.
<details>
<summary>Installation steps</summary>
<br>

```
sudo apt-get install cmake g++ make
```
</details>

### GSL
-- A library for numerical calculations, including matrix-vector multiplication, vector-vector operations, and random number generation.

<details>
<summary>Installation steps</summary>
<br>

Download the current stable version of the .tar.gz file from the gsl [webpage](https://www.gnu.org/software/gsl/).
Navigate to the directory where you downloaded the .tar.gz file, then:

1) extract the .tar.gz file:
```
tar zxvf gsl-latest.tar.gz
```
2) cd to the extracted directory
```
cd gsl-{version_number}/
```
3) check whether you have everything ready to build the application:
```
./configure
```
4) build the source code:
```
make
```
5) install it:
```
sudo make install
```
6) link the library - [docs](https://stackoverflow.com/questions/480764/linux-error-while-loading-shared-libraries-cannot-open-shared-object-file-no-s):
```
sudo ldconfig
```
</details>




## Installing this project
clone this project:
```
git clone <project_url>
```
cd inside the project folder:
```
cd PairPotMCinator/
```
create a build/ directory
```
mkdir build/
```
cd inside it
```
cd build/
```
run cmake to generate a makefile
```
cmake ../
```
build the source code using the makefile
```
make
```
install it
```
sudo make install
```

## Running the PairPotMCinator

Then cd to wherever you have your configuration .json file, and run
```
PairPotMCinator <your_json_file>
```

A few example .json files are located in the results/artificial/ folder of this project.

To run one of them, just cd there
```
cd results/artificial
```
and run it
```
PairPotMCinator config_simple.json
```

This will create a results directory `results_simple/` - storing the results.

For information regarding the structure of the json file and the storing of the results, refer to the `reference.md` file located in the root of this project.



## Editing the source code

a) Make your changes to the source code.

b) Recompile the source code:

cd to the build/ directory
```
cd build/
```
then run
```
cmake --build .
sudo make install
```