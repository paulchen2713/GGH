% 
% LLL original algorithm
%   ref. LLL algorithm pseudocode, wikipedia.
% 
clear;
clc;
%
% INPUT
%   a lattice basis B = {b0, b1, ... , bn} in Z^m
%   a parameter £_ with 0.25 < £_ < 1, most commonly £_ = 0.75 (3/4)
%
B = [19 2  32 46 3  33; 
     15 42 11 0  3  24; 
     43 15 0  24 4  16; 
     20 44 44 0  18 15; 
     0  48 35 16 31 31; 
     48 33 32 9  1  29];
%
% PROCEDURE
%
B_size = size(B, 1) - 1; % n == B_size
BB = GramSchmidt(B);     % compute B* and do not normalize
k = 1;
% while (k less or equal to n)
while k <= B_size
    % for (j form k-1 to 0)
    for j = k-1 : -1 : 0
        mu_kj = (B(k+1, :) * BB(j+1, :)') / (BB(j+1, :) * BB(j+1, :)'); % £g_i,j
        % if (|£g_k,j| > 1/2) then
        if abs(mu_kj) > 0.5
            B(k+1, :) = B(k+1, :) - round(mu_kj)*B(j+1, :);
            % update B* and the related £g_i,j 
            % (the naive method is to recompute B* whenever b_i changes)
            BB = GramSchmidt(B);
        end
    end
    mu_kk1 = (B(k+1, :) * BB(k, :)') / (BB(k, :) * BB(k, :)'); % £g_k,k-1
    bb_k = BB(k+1, :) * BB(k+1, :)'; % (b*_k)^2
    bb_k1 = BB(k, :) * BB(k, :)';    % (b*_k-1)^2
    % if (<b*_k, b*_k> greater or equal to (£_ - (£g_k,k-1)^2)<b*_k-1, b*_k-1> then
    if bb_k >= ((3/4) - mu_kk1^2) * bb_k1
        k = k + 1;
    else
        % swap b_k and b_k-1
        temp = B(k+1, :);
        B(k+1, :) = B(k, :);
        B(k, :) = temp;
        % updated B* and the related £g_i,j 
        % (the naive method is to recompute B* whenever b_i changes)
        BB = GramSchmidt(B);
        k = max(k-1, 1);
    end
end
% return B the LLL reduced basis of {b0, ... , bn}
%
% OUTPUT
%   the reduced basis b0, b1, ... , bn in Z^m
%
%
% testing through Example 7.75 M^LLL
%
M_LLL = [  7 -12  -8   4 19   9; 
         -20   4  -9  16 13  16;
           5   2  33   0 15  -9;
          -6  -7 -20 -21  8 -12;
         -10 -24  21 -15 -6 -11;
           7   4  -9 -11  1  31 ];
%
%
display(B);
% for in = 1 : B_size+1
%     if in == 1
%         fprintf('B = [ %d  %d  %d  %d  %d  %d;\n', B(in, :));
%     elseif in > 1 && in < B_size+1
%         fprintf('      %d  %d  %d  %d  %d  %d;\n', B(in, :));
%     elseif in == B_size+1
%         fprintf('      %d  %d  %d  %d  %d  %d ]\n\n', B(in, :));
%     end
% end
%
display(M_LLL);
% for in = 1 : B_size+1
%     if in == 1
%         fprintf('M_LLL = [ %d  %d  %d  %d  %d  %d;\n', M_LLL(in, :));
%     elseif in > 1 && in < B_size+1
%         fprintf('          %d  %d  %d  %d  %d  %d;\n', M_LLL(in, :));
%     elseif in == B_size+1
%         fprintf('          %d  %d  %d  %d  %d  %d ]\n\n', M_LLL(in, :));
%     end
% end
%
if isequal(B, M_LLL)
    fprintf('Correct, B == M_LLL\n\n');
else
    fprintf('Incorrect, B != M_LLL\n\n');
end
%
%

