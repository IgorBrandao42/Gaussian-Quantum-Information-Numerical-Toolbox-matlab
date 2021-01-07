classdef optical_cavity < handle         % class definning an optical cavity
  properties
    Delta                                % Natural frequency for the optical mode (Tweezer-cavity detunning in the Coherent Scattering setup)
    kappa                                % Cavity linewidth
    
    Entropy                              % Single mode entropy for cavity field
    nbar                                 % Occupation number   for cavity field
  end
  
  methods
    function obj = optical_cavity(Delta_0, kappa_0)
      % Class simulating an optical cavity interacting with N nanoparticles
      % The interaction is considered to be linear in the modes' positions quadratures
      % The cavity field is initially in a vacuum state
      %
      % PARAMETERS:
      %   Delta_0        - Natural frequency for the optical mode
      %   kappa_0        - Cavity linewidth
      
      obj.Delta = Delta_0;
      obj.kappa = kappa_0;
    end
    
  end
  
end

