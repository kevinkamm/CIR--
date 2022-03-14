function [Az,Bz] = RiccatiCIRTerminal(phi1z,phi2z,phi3z,az,bz,tau)
%%RICCATICIRTERMINAL computes the Riccati solutions of the CIR model for 
% given parameters
%    Input:
%       phi1z (double) : contains $\phi_1^z$
%       phi2z (double) : contains $\phi_2^z$
%       phi3z (double) : contains $\phi_3^z$
%       az (double): terminal value Az(T,T)=az
%       bz (double): terminal value Bz(T,T)=bz
%       tau (nd array): contains the time to maturities (T-t0)
%    Output:
%       Az (size(tau) array): contains the Riccati coefficient Az evaluated at
%                       T-t0, where t0 is the initial time and
%                       T the maturities 
%       Bz (size(tau) array): contains the Riccati coefficient Bz evaluated at
%                       T-t0, where t0 is the initial time and
%                       T the maturities 
%
temp=exp(phi1z.*tau)-1;
temp2=phi2z.*temp;
temp3=(1+bz.*(phi1z-phi2z));
Az = az.*(phi1z.*exp(phi2z.*tau)./...
     (phi1z + temp2.*temp3)).^phi3z;
Bz = (bz.*phi1z+temp.*temp3)./...
     (phi1z + temp2.*temp3);    
end