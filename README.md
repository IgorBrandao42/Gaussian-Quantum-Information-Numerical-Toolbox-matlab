# Numerical Linear Optomechanics

NuLO is a MATLAB Toolbox for time evolution of optomechanical systems

## Installation

Download this Toolbox and add its main folder to the MATLAB path:

```MATLAB
addpath('<download-path>/NuLO');
```

## Usage



#### Running Example
```c
addpath("<download-path>/NuLO)
Example_using_class
```

## Author
[Igor Brandão](mailto:igorbrandao@aluno.puc-rio.br) - Master's student in [Thiago Guerreiro](mailto:barbosa@puc-rio.br)'s Lab at Pontifícia Universidade Católica do Rio de Janeiro, Brazil

## Mathematical Formalism
For the optomechanics time evolution, this codes was originally created for and uses the same formalism as:
> Igor Brandão, Daniel Tandeitnik, Thiago Guerreiro, "Coherent Scattering-mediated heat transport between levitated nanospheres", [arXiv:2006.02455](https://arxiv.org/abs/2006.02455)

For the study of Gaussian quantum states, this code was based on and uses the same formalism as:

> Christian Weedbrook, Stefano Pirandola, Raúl García-Patrón, Nicolas J. Cerf, Timothy C. Ralph, Jeffrey H. , "Gaussian quantum information", [Rev. Mod. Phys. 84, 621](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.84.621)

## License
This code is made available under the Creative Commons Attribution - Non Commercial 4.0 License. For full details see LICENSE.md.

Cite this toolbox as: 
> Igor Brandão, "Numerical Linear Optomechanics" [https://github.com/IgorBrandao42/Linear-Quantum-Optomechanics](https://github.com/IgorBrandao42/Linear-Quantum-Optomechanics). Retrived <em>date you downloaded</em>

## File Listing
```C
simulation.m              - Class definning a simulation of N particles interacting through Coherent Scattering with a single cavity field mode 
particle.m                - Class definning a levitated nanoparticle
optical_cavity.m          - Class definning an optical cavity
Example_using_class.      - Example of usage
lyapunov_ode.m            - Numerical Lyapunov ODE for time evolution of the covariance matrix
logarithmic_negativity2.m - Calculation of the logarithmic negativity for a bipartite system
von_Neumann_Entropy2.m    - Calculation of the von Neumann entropy for a bipartite system
von_Neumann_Entropy1.m    - Calculation of the von Neumann entropy for a single mode
```
