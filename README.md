# DLBFoam: Dynamic load balancing for fast reactive simulations
![v1.1](https://img.shields.io/badge/DLBFoam-v1.1-blue)
![OpenFOAM 8](https://img.shields.io/badge/OpenFOAM-8-brightgreen)

## What is DLBFoam?
DLBFoam is an open-source library for OpenFOAM. It introduces dynamic load balancing and a zonal reference mapping model 
for fast chemistry calculation in parallel simulations. In addition, it also introduces a fully analytical Jacobian formulation and optimized ODE solution routines for further speed-up.

 
## Why do I need this?

Load imbalance in parallel reactive simulations is an issue that causes very long
simulation times in OpenFOAM simulations utilizing finite-rate chemistry.


![crab pet](https://i.imgur.com/yYVBgHV.gif)

## Prerequisites
- OpenFOAM installation (with correct version)

## Compilation

DLBFoam can be compiled by typing the following command after sourcing appropriate OpenFOAM version and making sure a valid LAPACK installation exists:

```
./Allwmake
```


## Usage

Once the compilation is successful, any case running with standard OpenFOAM can be easily converted to
use DLBFOAM, following these steps:

* The DLBFoam should be linked to the solver. Add the following to your system/controlDict file:

```
libs
(
    "libchemistryModel_DLB.so" 
);
```

* Select chemistry solver as ```ode_pyJac``` and the method as ```loadBalanced_pyJac``` in constant/chemistryProperties:

```
chemistryType
{
    solver          ode;
    method          loadBalanced;
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
I. Morev, B. Tekgül, M. Gadalla, A. Shahanaghi, J. Kannan, S. Karimkashi, O. Kaario, V. Vuorinen, Fast reactive flow simulations using analytical Jacobian and dynamic load balancing in OpenFOAM, arXiv preprint [arXiv:2105.12070](https://arxiv.org/abs/2105.12070) (2021).
<details>
<summary>BibTex</summary>
<p>
 
```
@article{morev2021fast,
  title={Fast reactive flow simulations using analytical Jacobian and dynamic load balancing in OpenFOAM},
  author={Morev, Ilya and Tekg{\"u}l, Bulut and Gadalla, Mahmoud and Shahanaghi, Ali and Kannan, Jeevananthan and Karimkashi, Shervin and Kaario, Ossi and Vuorinen, Ville},
  journal={arXiv preprint arXiv:2105.12070},
  year={2021}
}
```
 
</p>
</details>
