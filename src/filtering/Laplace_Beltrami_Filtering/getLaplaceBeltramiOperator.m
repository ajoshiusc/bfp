%
% LBO = getLaplaceBeltramiOperator(sf, precision)
% 
% Description:
%     Get Laplace Beltrami operator
% 
% Input:
%     sf - the surface structure
%     precision - precision of the operator
% 
% Output:
%     LBO - Laplace Beltrami operator structure
% 
% Copyright:
%     2016-2017 (c) USC Biomedical Imaging Group (BigLab)
% Author:
%     Anand Joshi, Jian (Andrew) Li
% Revision:
%     1.0.2
% Date:
%     2017/06/01
%

function LBO = getLaplaceBeltramiOperator(sf, precision)

    if ~exist('precision', 'var') || isempty(precision)
        precision = 0.5;
    end
    
    numV = size(sf.vertices, 1);
    
    if (precision > 0) && (precision <= 1)
        p = round(numV .* precision);
    elseif (precision >= 1) && (precision <= numV)
        p = round(precision);
    else
        error('incorrect precision');
    end
    
    SMat = get_stiffness_matrix_tri_wt(sf, ones(numV, 1));
    
    [E, Lambda] = eigs(SMat, p, 'sm');
    L = diag(Lambda);
    
    LBO.E = E;
    LBO.L = L;
    
end