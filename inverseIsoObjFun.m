% Vincente Pericoli
% UC Davis
% 4 Nov 2015
%
%Partially based on:
%   Silva, G. H., Le Riche, R., Molimard, J., & Vautrin, A. (2009). 
%   "Exact and efficient interpolation using finite elements shape 
%   functions." European Journal of Computational Mechanics/Revue 
%   Européenne de Mécanique Numérique, 18(3-4), 307-331.
%


function [ xi_, eflag ] = inverseIsoObjFun( point, nodCoord, LMTYPE )
%Inverse Isoparametric Mapping Objective Function Optimization
%using forward mapping and a Subspace Trust Region algorithm

%use a persistent variable for MATLAB optimization options, because
%setting these options takes a significant amount of time... so running
%this in a loop will take much longer if MYOPTIONS is reset each time
persistent MYOPTIONS

if isempty(MYOPTIONS)
    MYOPTIONS = optimoptions(@fminunc, 'Algorithm','trust-region', ...
                                       'MaxIter',1000, 'display','off', ...
                                       'GradObj','on', 'Hessian','on',  ...
                                       'TolX',eps(1), 'TolFun',eps(1));
end

[xi_, ~, eflag] = fminunc( @(xi_) ...
                    forwardMappingResid(xi_,point,nodCoord,LMTYPE), ...
                    [0,0,0], MYOPTIONS );


end


function [err, grad, Hess] = forwardMappingResid ...
                                          ( xi_, point, nodCoord, LMTYPE )
%forwardMappingResid - squared sum of the forward mapping residual

%obtain basis functions
[N, dN, d2N] = elementBasisFns( xi_, LMTYPE );

%calculate a squared "residual"
resid = [point(1) - N' * nodCoord(:,1);
         point(2) - N' * nodCoord(:,2);
         point(3) - N' * nodCoord(:,3)].^2;

%add up all points for a total "error"
err = resid(1) + resid(2) + resid(3);

%to use an optimization algorithm (subspace trust region), we also need the
%gradient and hessian of err

%err gradient:
grad = zeros(1,3);
for i = 1:3
    %for all derivative directions
    for p = 1:3
        %for all point coords x,y,z
        grad(i) = grad(i) ...
              - 2*(point(p) - N'*nodCoord(:,p))*(dN(:,i)' * nodCoord(:,p));
    end
end

%err Hessian:
%H(i,j) = f_,ij
Hess = zeros(3,3);
for i = 1:3
    %for all first derivative directions
    for j = 1:3
        %for all second derivative directions
        for p = 1:3
            %for all point coords x,y,z
            Hess(i,j) = Hess(i,j) + ...
              2*(dN(:,j)' * nodCoord(:,p))*(dN(:,i)' * nodCoord(:,p)) - ...
              2*(point(p) - N'*nodCoord(:,p))*(d2N(:,i,j)'*nodCoord(:,p));
        end
    end
end

return;
end
