function [phip, phiv] = calcPhipDistribution(C, C0, x)
% calulates linear solid fraction distribution across the pipe length
% Input: 
%     C - function slope [-]
%     C0 - constant (void fraction at the inlet, x=0) [-]
%     x - pipe length per particle dimater [-]
% Output: 
%    phip - solid fraction [-]
%    phiv - void fraction [-]


    % Element multiplication of a row vector and a column vector
    % results in a matrix. 
    [rows, columns] = size(x); % ensure result is a column vector
    if (rows == 1)
       x = x'; 
    endif

    phip = 1 - C0 - C .* x;  % solid fraction per length
    phiv = C0 + C .* x;      % void fraction per length

end

