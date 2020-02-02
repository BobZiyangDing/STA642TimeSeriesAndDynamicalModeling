function  [phi, ins]  = localstnytvar(p,phi,rho,A)                      
% Checks to see of an AR(p) coefficicent vector corresponds to stationarity
% If not, modify it (via eignevalues) to be stationary
% Inputs:
%   p    --  univariate TVAR model order
%   phi  --  px1 vector of TVAR coefficients  
%   rho  --  a number just less than 1 -- the max abs(eigenvalue); e.g., 0.99
%   A = [eye(p-1) zeros(p-1,1)];
% Outputs: 
%  phi   --  px1 vector of TVAR coefficients 
%  ins   -- 1 if the input phi was nonstationary, 0 otherwise 
 
ins=0; 
eigs = eig([ phi' ; A ]); ab = abs(eigs); ie=find(ab>=1);      % find if local non-stationarity? 
if (length(ie)>0)   % correct some eigenvalues
     ins=1; neweigs=eigs; ab(ie) = rho; 
     [x,y]=pol2cart(angle(eigs(ie)),ab(ie)); neweigs(ie)=complex(x,y); 
     %[ real(eigs) real(neweigs)  imag(eigs) imag(neweigs) ]
     %[ abs(eigs) abs(neweigs) angle(eigs)  angle(neweigs) ]
     phi = -real(poly(neweigs))'; phi(1)=[];    
end
    
    
    
