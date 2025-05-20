# PairPotMCinator

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
Then cd to wherever you have your .json, and run:
```
PairPotMCinator <your_json_file>
```

An example .json files are located in the results/artificial/ folder of this project.



## Hacking on the source code

a) Make your changes to the source code.

b) Recompile the source code:

cd to the build/ directory, then run
```
cmake --build .
sudo make install
```