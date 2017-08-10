%
% [x, output] = fastNNLS(AtA, Atb, optNNLS)
% 
% Description:
%     fast non-negative least square implementation of Bro's paper
% 
% Input:
%     AtA - A' * A
%     Atb - A' * b
%     optNNLS - options specify tol and max iteration
% 
% Output:
%     x - solution
%     output - tol and actual number of iterations
% 
% Copyright:
%     2016-2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Jian (Andrew) Li
% Revision:
%     2.3.2
% Date:
%     2017/06/01
%

function [x, output] = fastNNLS(AtA, Atb, optNNLS)

    if nargin == 0 % without any parameter return default option
        optNNLS = struct;
        optNNLS.tol = [];
        optNNLS.maxNumItrOuter = 1e50;
        optNNLS.maxNumItrInner = [];
        x = optNNLS;
        output = [];
        return;
    end
    
    if ~exist('option', 'var') || isempty(optNNLS)
        optNNLS = fastNNLS();
    end
    
    if isempty(optNNLS.tol)
        tol = 10 * eps * norm(AtA, 1) * max(size(AtA));
    else
        tol = optNNLS.tol;
    end
    
    N = size(AtA, 2);
    wz = zeros(N, 1);

    P = false(N, 1);
    Z = true(N,1);
    x = zeros(N, 1);

    w = Atb - AtA*x;

    cIn = 0;
    cOut = 0;
    
    if isempty(optNNLS.maxNumItrInner)
        maxNumItrInner = 3 * N;
    else
        maxNumItrInner = optNNLS.maxNumItrInner;
    end
    
    maxNumItrOuter = optNNLS.maxNumItrOuter;

    while any(Z) && any(w(Z) > tol) && (cOut < maxNumItrOuter)
        cOut = cOut + 1;
        wz(P) = -Inf;
        wz(Z) = w(Z);

        [~, t] = max(wz);
        P(t) = true;
        Z(t) = false;

        z = zeros(N, 1);
        z(P) = AtA(P, P) \ Atb(P);

        while any(z(P) <= 0) && (cIn < maxNumItrInner)
%             disp('in');
            cIn = cIn + 1;
            Q = (z <= 0) & P;
            alpha = min(x(Q)./(x(Q) - z(Q)));
            x = x + alpha*(z - x);
            Z = ((abs(x) < tol) & P) | Z;
            P = ~Z;
            z = zeros(N, 1);
            z(P) = AtA(P, P) \ Atb(P);
        end
        x = z;
        w = Atb - AtA*x;
    end
    
    output = struct;
    output.tol = tol;
    output.numItr = cOut;  
end