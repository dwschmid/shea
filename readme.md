# shea
Shear HEAting model for 1d shear zone implemented in Matlab (v2020a).

**shea** was developed as a contribution to "Rapid fluid-driven transformation of lower continental crust associated with thrust-induced shear heating" by Jamtveit et al. (2021), published in Lithos [doi:10.1016/j.lithos.2021.106216](https://www.sciencedirect.com/science/article/pii/S0024493721002528?dgcid=coauthor). 

* **shea_compute.m**      - main computational routine
* **shea_constants.m**    - definition of various constants
* **shea_driver.m**       - driver for shea_compute, facilitates parameter space analysis
* **shea_materials.m**    - database of viscous flow law and thermal parameters 
* **shearzone_nf.m**      - shear zone model with no thermo-mechanical feedback

Archived on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3994114.svg)](https://doi.org/10.5281/zenodo.3994114)

# Background
## Governing Equations
Our implementation is a coupled thermal and viscous mechanical model in 1D. 

We assume that the model is aligned with spatial coordinate ![z](https://render.githubusercontent.com/render/math?math=z) and that the shearing takes place parallel to ![x](https://render.githubusercontent.com/render/math?math=x). The governing equations are linear momentum:

![\frac{{\partial {\tau _{xz}}}}{{\partial z}} = 0](https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%7B%5Cpartial%20%7B%5Ctau%20_%7Bxz%7D%7D%7D%7D%7B%7B%5Cpartial%20z%7D%7D%20%3D%200)

and energy conservation

![\rho {c_p}\frac{{\partial T}}{{\partial t}} = \frac{\partial }{{\partial z}}\left( {\lambda \frac{{\partial T}}{{\partial z}}} \right) + {\tau _{xz}}{\dot \varepsilon _{xz}}](https://render.githubusercontent.com/render/math?math=%5Crho%20%7Bc_p%7D%5Cfrac%7B%7B%5Cpartial%20T%7D%7D%7B%7B%5Cpartial%20t%7D%7D%20%3D%20%5Cfrac%7B%5Cpartial%20%7D%7B%7B%5Cpartial%20z%7D%7D%5Cleft(%20%7B%5Clambda%20%5Cfrac%7B%7B%5Cpartial%20T%7D%7D%7B%7B%5Cpartial%20z%7D%7D%7D%20%5Cright)%20%2B%20%7B%5Ctau%20_%7Bxz%7D%7D%7B%5Cdot%20%5Cvarepsilon%20_%7Bxz%7D%7D)

Here ![${\tau _{xz}}$](https://render.githubusercontent.com/render/math?math=%24%7B%5Ctau%20_%7Bxz%7D%7D%24) is the shear stress, ![\rho](https://render.githubusercontent.com/render/math?math=%5Crho) density, ![\c_p](https://render.githubusercontent.com/render/math?math=%5Cc_p) heat capacity,  ![T](https://render.githubusercontent.com/render/math?math=T) temperature, ![t](https://render.githubusercontent.com/render/math?math=t) time, and ![\lambda](https://render.githubusercontent.com/render/math?math=%5Clambda) the conductivity. 
The shear strain rate ![\dot \varepsilon _{xz}](https://render.githubusercontent.com/render/math?math=%5Cdot%20%5Cvarepsilon%20_%7Bxz%7D) is related to the shear velocity ![v_x](https://render.githubusercontent.com/render/math?math=v_x) by

![{\dot \varepsilon _{xz}} = \frac{1}{2}\frac{{\partial {v_x}}}{{\partial z}}](https://render.githubusercontent.com/render/math?math=%7B%5Cdot%20%5Cvarepsilon%20_%7Bxz%7D%7D%20%3D%20%5Cfrac%7B1%7D%7B2%7D%5Cfrac%7B%7B%5Cpartial%20%7Bv_x%7D%7D%7D%7B%7B%5Cpartial%20z%7D%7D)

The constitutive equation is 

![{\tau _{xz}} = 2\eta {\dot \varepsilon _{xz}}](https://render.githubusercontent.com/render/math?math=%7B%5Ctau%20_%7Bxz%7D%7D%20%3D%202%5Ceta%20%7B%5Cdot%20%5Cvarepsilon%20_%7Bxz%7D%7D)

where ![\eta](https://render.githubusercontent.com/render/math?math=%5Ceta) is the effective viscosity:

![\eta  = F{A^{ - \frac{1}{n}}}\dot \varepsilon _{II}^{\frac{1}{n} - 1}{e^{\frac{Q}{{nRT}}}}](https://render.githubusercontent.com/render/math?math=%5Ceta%20%20%3D%20F%7BA%5E%7B%20-%20%5Cfrac%7B1%7D%7Bn%7D%7D%7D%5Cdot%20%5Cvarepsilon%20_%7BII%7D%5E%7B%5Cfrac%7B1%7D%7Bn%7D%20-%201%7D%7Be%5E%7B%5Cfrac%7BQ%7D%7B%7BnRT%7D%7D%7D%7D)

Here ![F](https://render.githubusercontent.com/render/math?math=F) is a geometrical conversion factor, ![A](https://render.githubusercontent.com/render/math?math=A) is the pre-exponential factor of the flaw law, ![Q](https://render.githubusercontent.com/render/math?math=Q) the activation energy, ![n](https://render.githubusercontent.com/render/math?math=n) the power law exponent, and ![R](https://render.githubusercontent.com/render/math?math=R) the universal gas constant. ![{\dot \varepsilon _{II}}](https://render.githubusercontent.com/render/math?math=%7B%5Cdot%20%5Cvarepsilon%20_%7BII%7D%7D) is the is the square root of the second invariant of the strain rate tensor which in our case simplifies to

![{\dot \varepsilon _{II}} = {\dot \varepsilon _{xz}}](https://render.githubusercontent.com/render/math?math=%7B%5Cdot%20%5Cvarepsilon%20_%7BII%7D%7D%20%3D%20%7B%5Cdot%20%5Cvarepsilon%20_%7Bxz%7D%7D)

## Geometry and Boundary Conditions
We specify the depth to the top of the shear zone and the initial temperature in the center of the the shear zone. 
These values are based on independent pressure and temperature estimates. 
We then assume that the surface boundary condition temperature is 0 and based on this we initialize the intitial temperature with a constant gradient. 
The temperature model comprises the overburden, the shear zone, and host rock under the shear zone with a thickness equivalent to the overburden. 
The bottom boundary condition is set to match the initial temperature gradient. 
The material properties for the entire model are set to represent a single lithology and we ignore the potential heat source contribution due to radioactive decay. 
The resulting initial geotherm does not aim to represent a realistic crustal geotherm but is designed to match the inferred background state of the shear zone prior to shearing
and the size of the model is kept large as to reduce the influence of effect of temperature boundary conditions.

The mechanical model only solves for the deformation in the shear zone. 
The boundary condition is a constant displacement velocity. 
For a given temperature state the mechanical solver iteratively determines the velocity profile throughout the shear zone and the related effective viscosity. 
The energy dissipation due to the shearing is then accounted for in the next temperature solve. 

## Numerical Method
We numerical method is a finite difference method. The mechnical part is solved implicitely while we use an explicite solver for the temperature. 

## Examples
Running **shea_compute** shows how the various parameters (temperature, viscosity, etc.) evolve for a specific configuration of the controlling parameters. The result visualization is also saved as an animated gif, which we show below:
![Parameter evolution in shear zone](/doc/shea_example.gif?raw=true)

Running **shea_driver.m** scans through a range of parameters and records the maximum stress and temperature values experienced in the center of the shear zone. An example is shown below:
![Maximum temperature in shear zone for parameter ranges](/doc/parameter_study_temp.png?raw=true)

## Acknowledgements
This project was supported by the European Union's Horizon 2020 Research and Innovation Programme under the ERC Advanced Grant Agreement n°669972, ‘Disequilibrium Metamorphism’ (‘DIME’) to BJ.
