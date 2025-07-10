function Pdf=MGGD_fun(X,Sigma,Beta,mu,m,p)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2024
%   DESCRIPTION:

%   This function gets an observed data vector and a set of MGGD parameters
%   and outputs its propability density function (pdf).

%   INPUT:
%   > X                   - 3XM Matrix containing the 3D features of M observed data points
%   > Sigma               - 3X3 Covariance matrix 
%   > Beta                - 3X1 Shape vector 
%   > mu                  - 3X1 Expectaion Vector 
%   > p                   - 3X1 Scale vector 

%   OUTPUT:
%   > Pdf                 - Vector of MX1 representing the pdf of all input data points.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=(gamma(p/2))/((pi^(p/2))*gamma(p/(2*Beta))*(2^(p/(2*Beta))));
B=Beta/((m^(p/2))*(sqrt(det(Sigma))));
a=1/(2*(m^Beta));
if p==3
    MU=[mu(1)*ones(size(X,1),1) mu(2)*ones(size(X,1),1) mu(3)*ones(size(X,1),1)];
elseif p==2
    MU=[mu(1)*ones(size(X,1),1) mu(2)*ones(size(X,1),1)];
end
b=(X-MU)*inv(Sigma)*(X-MU)';
C=exp(-a*(b.^Beta));

Pdf=A*B*C;

end










