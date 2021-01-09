
# Gaussian Quantum Information Toolbox for Linear Optomechanics

This is a MATLAB Toolbox for numerical simulation of the time evolution of linear quantum optomechanical systems in gaussian states following an open quantum dynamics.

The simulated system consists of N particles interacting with a single eletromagnetic field, considered to be an optical cavity field. The Hamiltonian that describes this system is

![equation](https://latex.codecogs.com/gif.latex?%5Chat%7BH%7D%20%3D%20%5Chbar%5CDelta%20%5Chat%7Ba%7D%5E%5Cdagger%5Chat%7Ba%7D%20&plus;%20%5Csum_%7Bj%3D1%7D%5E%7BN%7D%20%5Chbar%5Comega_j%20%5Chat%7Bb%7D_j%5E%5Cdagger%5Chat%7Bb%7D_j%20&plus;%20%5Csum_%7Bj%3D1%7D%5E%7BN%7D%20%5Chbar%20g_j%28%5Chat%7Ba%7D%5E%5Cdagger%20&plus;%20%5Chat%7Ba%7D%29%28%5Chat%7Bb%7D_j%5E%5Cdagger%20&plus;%20%5Chat%7Bb%7D_j%29%20%5C%2C%20%2C)

where ![equation](https://latex.codecogs.com/gif.latex?%5Chat%7Ba%7D) (![equation](https://latex.codecogs.com/gif.latex?%5Chat%7Ba%7D%5E%5Cdagger)) is the annihilation (creation) operator for the cavity field with natural frequency ![equation](https://latex.codecogs.com/gif.latex?%5CDelta); ![equation](https://latex.codecogs.com/gif.latex?%5Chat%7Bb%7D_j) (![equation](https://latex.codecogs.com/gif.latex?%5Chat%7Bb%7D%5E%5Cdagger_j)) is the annihilation (creation) operator for the j-th particle each with natural frequency ![equation](https://latex.codecogs.com/gif.latex?%5Comega_j), and ![equation](https://latex.codecogs.com/gif.latex?g_j) is the coupling strength of the interaction.

Initially, the particles are assumed to be each in a thermal state and the cavity field in a vacuum state. As this Hamiltonian preserves their gaussianity [[Rev. Mod. Phys. 84, 621]](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.84.621), this toolbox focus on the time evolution of the expected values of the quadratures and the covariance matrix.

The time evolution is modelled by a set of quantum Langevin equations for the quadratures and a Lyapunov equation for the covariance matrix. For the full description of the open quantum dynamics, please refer to Section III of [[Brandao2021]](https://arxiv.org/abs/2006.02455).

## Installation

Clone this repository or download this Toolbox and add its main folder to the MATLAB path:

```MATLAB
addpath('<download-path>/<name-folder>');
```

## Usage

This program only inputs are the parameters values and time interval for the calculation. The program infers the number of particles from the lengths of these vectors.
```MATLAB
% Define the constants for the system under study. 
% In this case 3 particles and one cavity.
omega =    2*pi*[190e+3;  160e+3;  180e+3];     % Natural frequency of the particles
g     =    2*pi*[ 42e+3; 35.3e+3; 39.8e+3];     % Coupling strength
gamma =    2*pi*[ 10e+3;   10e+3;   10e+3];     % Damping
T     =         [  1e-1;    1e-1;    1e-1];     % Initial temperature of each particle
T_environment = [  1e-1;    1e-1;    1e-1];     % Temperature for the environment of each particle
Delta = +omega(1);                              % Cavity field natural frequency
kappa = 2*pi*193e+3;                            % Cavity linewidth

% Define the time interval to be studied
t = linspace(0, 1e-5, 1e+3);
```

You need to create an instance of a simulation (class) with the parameters and run the calculations at the time stamps:
```MATLAB
% Create an instance of a simulation
example = simulation(omega, g, gamma, T, T_environment, Delta, kappa);

% Run the every calculation available
example.run(t);
```

The simulations are a handle class and its results can be directly extracted or plotted using its internal plotting methods:
```MATLAB
% Plot the results
example.plot();
```

The user can choose to calculate only what suits them, by passing extra parameters to the method 'run', they are:

|  Optional parameter  | Specific calculation of time evolution |
|----------------------|----------------------|
| "langevin"           | Solve semiclassical Langevin equations for the expectation value of the quadratures |
| "lyapunov"           | Solve Lyapunov equation for the covariance matrix| 
| "steady_state"       | Find steady state covariance matrix |
|"occupation_number"   | Find the occupation number for each mode|
| "entanglement"       | Calculate the logarithmic negativity for each bipartition |
| "entropy"            | Calculate the von Neumann entropy for each mode, bipartition and the whole system|
| "mutual_information" | Calculate the mutual information for the whole system|
|   "heat_flux"        | Calculate the expectation value of every heat fluxes in the system|
| "fidelity_test"      | Approximate each mode state by a thermal state through Fideliy, finding the effective temperature |

```MATLAB
% Run a specific calculation
example.run(t, "occupation_number");
```

#### Running Example
In the file **Example.m** there is a basic example of the capabilities of this Toolbox

## Author
[Igor Brandão](mailto:igorbrandao@aluno.puc-rio.br) - Master's student in [Thiago Guerreiro](mailto:barbosa@puc-rio.br)'s Lab at Pontifical Catholic University of Rio de Janeiro, Brazil

## Mathematical Formalism
For the optomechanics time evolution, this codes was originally created for and uses the same formalism as:
> Igor Brandão, Daniel Tandeitnik, Thiago Guerreiro, "Coherent Scattering-mediated heat transport between levitated nanospheres", [arXiv:2006.02455](https://arxiv.org/abs/2006.02455)

For the study of Gaussian Quantum Information, this code was based on and uses the same formalism as:

> Christian Weedbrook, Stefano Pirandola, Raúl García-Patrón, Nicolas J. Cerf, Timothy C. Ralph, Jeffrey H. , "Gaussian quantum information", [Rev. Mod. Phys. 84, 621](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.84.621)

## License
This code is made available under the Creative Commons Attribution - Non Commercial 4.0 License. For full details see LICENSE.md.

Cite this toolbox as: 
> Igor Brandão, "Gaussian Quantum Information Toolbox for Linear Optomechanics", [https://github.com/IgorBrandao42/Gaussian-Quantum-Information-Toolbox-for-Linear-Optomechanics](https://github.com/IgorBrandao42/Gaussian-Quantum-Information-Toolbox-for-Linear-Optomechanics). Retrived <*date you downloaded*>

## File Listing

|    File name     |              What it does |
|------------------|----------------------------------------------|
|  simulation.m    |  Class definning a simulation of N particles |
|   particle.m  |    Class definning a nanoparticle  |
|   optical_cavity.m  |   Class definning an optical cavity   |
|   symplectic_eigenvalues.m  |  Calculates the sympletic eigenvalues of a covariance matrix    |
|   von_Neumann_Entropy.m  |  Calculates the von Neumann entropy for a multipartite gaussian system from its covariance matrix  |
|   logarithmic_negativity2.m  |   Calculation of the logarithmic negativity for a bipartite system from its covariance matrix   |
|    single_mode_CM.m   |     Finds the covariance submatrix for a single mode from the full covariance matrix      |
|    bipartite_CM.m   |       Finds the covariance submatrix for a bipartition from the full covariance matrix     |
|    lyapunov_ode    | ODE that defines the Lyapunov equation to be integrated by ode45 |
  |  func.m    | Auxiliar mathematical function to the von Neumann entropy |
  | Example.m     | Basic example of usage of the Toolbox

## Acknowledgment
The author thanks Daniel Ribas Tandeitnik and Professor Thiago Guerreiro for the discussions. The author is thankful for support received from FAPERJ Scholarship No. E-26/200.270/2020



