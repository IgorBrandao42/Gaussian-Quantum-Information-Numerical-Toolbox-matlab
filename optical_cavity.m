classdef optical_cavity < handle                           % class definning an optical cavity
  properties
    Delta                          % Tweezer-cavity detunning
    kappa                          % Caity linewidth
    
    Entropy                        % Single mode entropy for cavity field
    nbar                           % Occupation number   for cavity field
  end
  
  methods
    function obj = optical_cavity(Delta_0, kappa_0)
      obj.Delta = Delta_0;
      obj.kappa = kappa_0;
    end
    
  end
  
end

