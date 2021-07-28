%% Creation of many different state
vacuum = gaussian_state("vacuum");

thermal = gaussian_state("thermal", 100);

squeezed = gaussian_state("squeezed", 1.2);
squeezed.displace(2 + 5j);

coherent = gaussian_state("coherent", 2+1j);
coherent.rotate(pi/6);

bipartite = squeezed.tensor_product([thermal]);
bipartite.two_mode_squeezing(1.5);

partition = bipartite.partial_trace([1]);

tripartite = thermal.tensor_product([squeezed, coherent]);

bipartition = tripartite.only_modes([1,2]);


%% Retrieval of information from the gaussian states
lambda_thermal = thermal.symplectic_eigenvalues();

lambda_tri = tripartite.symplectic_eigenvalues();

p_tri = tripartite.purity();

p_coherent = coherent.purity();

S = thermal.von_Neumann_Entropy();

I = tripartite.mutual_information();

nbar_th = thermal.occupation_number();

nbar_3 = tripartite.occupation_number();

F_ac = vacuum.fidelity(coherent);

C_squeezed = squeezed.coherence();


