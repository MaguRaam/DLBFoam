# DLBFoam: Dynamic load balancing for fast reactive simulations
![v1.1](https://img.shields.io/badge/DLBFoam-v1.1-blue)
![OpenFOAM dev](https://img.shields.io/badge/OpenFOAM-dev_20230829-brightgreen)

[<img src="https://img.shields.io/badge/-YouTube_Video_Tutorials-red?style=for-the-badge&logo=youtube&logoColor=white"/>](https://www.youtube.com/playlist?list=PLXqVaOXSsv1SBnfyGRa_C-E0X--FIT27P)

## DLBFoam v1.1 - What's new?
DLBFoam v1.1 introduces a fully analytical chemistry Jacobian via [pyJac](https://github.com/SLACKHA/pyJac), and optimized ODE solution routines via [LAPACK](http://www.netlib.org/lapack/). Combined with the load balancing features, v1.1 provides up to x250 speed-up compared to standard OpenFOAM chemistry model. If you are interested with using only dynamic load balancing without any third party dependency, please use [DLBFoam v1.0](https://github.com/blttkgl/DLBFoam-1.0).

## What is DLBFoam?
DLBFoam is an open-source library for OpenFOAM. It introduces dynamic load balancing and a zonal reference mapping model 
for fast chemistry calculation in parallel simulations. In addition, it also introduces a fully analytical Jacobian formulation and optimized ODE solution routines for further speed-up.

 
## Why do I need this?

Load imbalance in parallel reactive simulations is an issue that causes very long
simulation times in OpenFOAM simulations utilizing finite-rate chemistry.

DLBFoam introduces runtime load balancing through MPI routines
to minimize the load imbalance between ranks and gain speed-up. The implementation
details can be found in our paper [[1]](#1). In addition, the cell-wise chemistry problem is vastly improved by the analytical Jacobian formulation and optimized matrix operations in the ODE solver class. The details for those implementations can be found in our follow-up paper [[2]](#2).


![crab pet](https://i.imgur.com/yYVBgHV.gif)

## Prerequisites
- OpenFOAM installation (with correct version)
- LAPACK (Intel-MKL, OpenBLAS or standalone)
- Cmake
- [ct2foam](https://github.com/kahilah/ct2foam) (Optional)

## Compilation

DLBFoam can be compiled by typing the following command after sourcing appropriate OpenFOAM version and making sure a valid LAPACK installation exists:

```
./Allwmake --clean --platform <LAPACK_INSTALLATION_TYPE>
```
<LAPACK_INSTALLATION_TYPE> can be MKL, OPENBLAS or STANDALONE.


DLBFoam requires LAPACK packages for improved ODE routines (LAPACKE C interface for OPENBLAS and standalone installation). There are three different ways to provide LAPACK for DLBFoam:

- **Intel-MKL**:  We recommend Intel-MKL libraries to be used together with DLBFoam whenever you are working on machines with Intel-based architecture. See further information and installation guidelines for  [Intel-MKL](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html). Note that DLBFoam assumes, that the ```MKLROOT``` environment variable is set to represent the installation path according to the standard library installation scripts.

- **OpenBLAS**: Another option is to utilise [OpenBLAS](https://www.openblas.net/) library which includes LAPACK routines. In this case, DLBFoam assumes that the ```OPENBLAS_INSTALL_ROOT``` environment variable is set to represent the OpenBLAS installation path succesfully.

- **Standalone**: A standalone installation may be a good idea if you are on your personal workstation and not on a cluster. You can see if you have the necessary lapacke dependency by e.g. executing ```ldconfig -p | grep "liblapacke"```. In case not, on debian systems you could e.g. install the requirements (including header files) by:

    ```
    (sudo) apt-get install liblapacke-dev
    ```

## Mechanism generation
After the successful compilation, mechanism files for tutorials are compiled automatically, so you can run tutorials right away.
If you want to use another mechanism, the C subroutines generated by the pyJac for the analytical Jacobian should be compiled and linked to the OpenFOAM case folder. 

In this repository, we provide the input files and C subroutines generated by pyJac for [GRI-3.0](http://combustion.berkeley.edu/gri-mech/version30/text30.html), [Yao](https://www.sciencedirect.com/science/article/abs/pii/S001623611631184X), and [DRM-19](http://combustion.berkeley.edu/drm/) mechanisms. If you want to use one of these mechanisms, no other dependency is required (except a ```CMake``` installation).

* Go to the ```utilities/``` folder, and execute ```Allrun``` as:

```
./Allrun -m <mechName>
```
where ```<mechName>``` is ```gri30```, ```yao```, or ```drm19```. This command will generate a compiled mechanism library ```<mechName>/lib/build/libc_pyjac.so```. You need to put this mechanism library along with properties file located in ```<mechName>/foam/``` to e.g. ```constant/foam/``` folder and link them in ```controlDict```, as explained in the next section.

**If you want to use a different chemical mechanism**, you need to create the thermo input files in a format required by OpenFOAM, as well as analytical Jacobian C subroutines generated by pyJac. This process requires many dependencies, most notably [Cantera](https://cantera.org/) and pyJac.

We have developed a tool called
[ct2foam](https://github.com/kahilah/ct2foam) that makes this process a lot easier. Please check ct2foam (and the pyjac2foam utility in it) if you are interested in using a different mechanism.
## Usage

Once the compilation is successful, any case running with standard OpenFOAM can be easily converted to
use DLBFOAM, following these steps:

* The DLBFoam should be linked to the solver. Add the following to your system/controlDict file:

```
libs
(
    "libchemistryModel_DLB.so" 
    "libODE_DLB.so"
    "$FOAM_CASE/constant/foam/libc_pyjac.so"
);
```
the first two libraries link the DLBFoam and the optimized LAPACK solvers, while the last library links the C subroutines generated in the previous section for the analytical Jacobian.

* Select chemistry solver as ```ode_pyJac``` and the method as ```loadBalanced_pyJac``` in constant/chemistryProperties:

```
chemistryType
{
    solver          ode_pyJac;
    method          loadBalanced_pyJac;
}
```

* Add the loadbalancing subdictionary to the same chemistryProperties file:

```
loadbalancing
{
    active true;
    log	true;
}
```

* Set the solver flag under ```odeCoeffs``` to ```seulex_LAPACK``` in order to use the optimized ODE solvers:
```
odeCoeffs
{
    solver          seulex_LAPACK;
    absTol          1e-08;
    relTol          1e-05;
}
```
* (Optional) Set the refmapping as active in chemistryProperties file if you want to 
    use the reference mapping method (you have to add an empty ```refmapping{}``` dict
    even if you do not use it):

```
refmapping
{
    active  true;
    
    mixtureFractionProperties
    {
        oxidizerMassFractions
        {
            N2       0.77;
            O2       0.23;
        }

        fuelMassFractions
        {
            NC12H26       1.0;
        }

        #include "$FOAM_CASE/constant/foam/thermo.foam"
    }
    tolerance	1e-4;  // mixture fraction tolerance
    deltaT	2; // temperature tolerance
}
```
Reference mapping uses mixture fraction (Z) and maps a reference solution to reference
cells satisfying a condition.

The entry above sets the Z=0 and Z=1 conditions from given mass fractions. For each
CFD iteration it finds a reference solution where Z<tolerance and solves the chemistry.
Subsequent cells following the same condition are mapped from this reference solution.

(Optional) When deltaT is explicitly set, the mapper also checks the temperature
between reference solution and other reference cells and ensures:
abs(T<sub>cell</sub>-T<sub>ref</sub>)<deltaT.


* Run the case normally with OpenFOAM's reactive solvers.

For a working example, check the tutorials given in tutorials folder.

## FAQ

#### Undefined symbol: eval_h
```
/path/to/user/OpenFoam/linux64GccDPInt32Opt/lib/libchemistryModel_DLB.so: undefined symbol: eval_h
```
Compiled mechanism library (```libc_pyjac.so```) is not found. Please read Mechanism Generation section of the README. Check if a correct path is set in controlDict and the file exists.

#### Tutorial works in serial, but it hangs in parallel
The reason might be that both OpenBLAS and LAPACKE are installed in your system and interfering with each other. This mostly happens in personal computers rather than clusters where software are controlled by modules. On your personal system, please follow the instructions below to mitigate this issue:

1. Find all the OpenBlas related packages installed on your machine with ```dpkg --get-selections | grep openblas``` and remove all of them.
2. Remove ```liblapacke-dev``` and install it again as ```(sudo) apt-get install liblapacke-dev```.
3. Recompile the code with the compilation flag ```--clean``` for ```STANDALONE``` platform: ```./Allwmake --clean --platform STANDALONE```

## Contributors
- Bulut Tekgül (buluttekgul@gmail.com)
- Petteri Peltonen (petteri.peltonen@aalto.fi)
- Heikki Kahila (heikki.kahila@wartsila.com)
- Ilya Morev (ilya.morev@aalto.fi)
- Mahmoud Gadalla (mahmoud.gadalla@aalto.fi)


## Getting help and reporting bugs

Please submit a GitHub issue if you found a bug in the program. If you need help with the software or have further questions, either open an issue or contact the contributors.

## Citation

If you use our model, please cite publications describing its implementation, Refs. [[1]](#1) and [[2]](#2). 

## References

<a id="1">[1]</a> 
B. Tekgül,  P. Peltonen,  H. Kahila,  O. Kaario,  V. Vuorinen,  DLBFoam: An open-source dynamic load balancing model for fast reacting flow simulations in OpenFOAM, Computer Physics Communications, Volume 267, [10.1016/j.cpc.2021.108073](https://doi.org/10.1016/j.cpc.2021.108073) (2021).
<details>
<summary>BibTex</summary>
<p>
 
```
@article{tekgul2021dlbfoam,
  title={DLBFoam: An open-source dynamic load balancing model for fast reacting flow simulations in OpenFOAM},
  author={Tekg{\"u}l, Bulut and Peltonen, Petteri and Kahila, Heikki and Kaario, Ossi and Vuorinen, Ville},
  journal={Computer Physics Communications},
  pages={108073},
  year={2021},
  publisher={Elsevier}
}
```
 
</p>
</details>

<a id="2">[2]</a> 
I. Morev, B. Tekgül, M. Gadalla, A. Shahanaghi, J. Kannan, S. Karimkashi, O. Kaario, V. Vuorinen, Fast reactive flow simulations using analytical Jacobian and dynamic load balancing in OpenFOAM, Physics of Fluids 34, 021801, [10.1063/5.0077437](https://doi.org/10.1063/5.0077437) (2022).
<details>
<summary>BibTex</summary>
<p>
 
```
@article{morev2022fast,
  author = {Morev,Ilya  and Tekg{\"u}l,Bulut  and Gadalla,Mahmoud  and Shahanaghi,Ali  and Kannan,Jeevananthan  and Karimkashi,Shervin  and Kaario,Ossi  and Vuorinen,Ville },
  title = {{Fast reactive flow simulations using analytical Jacobian and dynamic load balancing in OpenFOAM}},
  journal = {Physics of Fluids},
  volume = {34},
  number = {2},
  pages = {021801},
  year = {2022},
  doi = {10.1063/5.0077437},
}
```
 
</p>
</details>
