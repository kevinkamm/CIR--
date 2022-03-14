function [Az,Bz] = RiccatiCIR(phi1z,phi2z,phi3z,tau)
%RICCATICIR computes the Riccati solutions of the CIR model for given parameters
%    Input:
%       phi1z (double) : contains $\phi_1^z$
%       phi2z (double) : contains $\phi_2^z$
%       phi3z (double) : contains $\phi_3^z$
%       tau (nx1 array): contains the time to maturities (T-t0)
%    Output:
%       Az (nx1 array): contains the Riccati coefficient Az evaluated at
%                       T-t0, where t0 is the initial time and
%                       T the maturities 
%       Bz (nx1 array): contains the Riccati coefficient Bz evaluated at
%                       T-t0, where t0 is the initial time and
%                       T the maturities 
%
temp=exp(phi1z.*tau)-1;
temp2=phi2z.*temp;
Az = (phi1z.*exp(phi2z.*tau)./...
     (phi1z + temp2)).^phi3z;
Bz = temp./...
     (phi1z + temp2);    
end