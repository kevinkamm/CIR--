function [A,varargout] = linconstr(varargin)
%LINCONSTR provides nonlinear constraints for CIR- model
%    Input:
%       varargin (cell array): contains (name,value) pairs
%           'b': its value contains the right-hand side of the 
%                linear inequality constraints
%                default: [0 0 0 0]
%    Output:
%       A (4x8 array): contains the linear constraints for given
%                      parameters in the form A*x<=b
%    Usage:
%       A=linconstr(varargin)
%       [A,b]=linconstr(varargin)
%
% See also calibrationGC, calibrationMC, fmincon

% b=[0,0,0,0,0.01,0.01];
b=[0,0,0,0];
for k=1:2:nargin
    switch varargin{k}
        case 'b'
            b = varargin{k+1};
    end
end

% % parameters for $x_t$
% phi1x = params(1);
% phi2x = params(2);
% phi3x = params(3);
% xt0   = params(7);
% 
% % parameters for $y_t$
% phi1y = params(4);
% phi2y = params(5);
% phi3y = params(6);
% yt0   = params(8);

% A=zeros(6,8);
A=zeros(4,8);

% % $\sigma_x \in \mathbb{R} \Leftrightarrow \phi_1^x \geq \phi_2^x$
% c(1) = phi2x - phi1x;
A(1,1:2)=[-1 1];

% % $\sigma_y \in \mathbb{R} \Leftrightarrow \phi_2^y \geq \phi_1^y$
% c(2) = phi1y - phi2y;
A(2,4:5)=[1 -1];

% % Positive mean-reversion speed $2\phi_2^x \geq \phi_1^x$
% c(3) = phi1x - 2*phi2x;
A(3,1:2)=[1 -2];

% % Positive mean-reversion speed $2\phi_2^y \geq \phi_1^y$
% c(4) = phi1y - 2*phi2y;
A(4,4:5)=[1 -2];

% A(5,7:8)=[1 -1];
% A(6,7:8)=[-1 1];

switch nargout
    case 2
        varargout{1} = b;
end
end