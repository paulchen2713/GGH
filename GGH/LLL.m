% 
% LLL algorithm
%   ref. LLL algorithm pseudocode, wikipedia.
% 
function B = LLL(B)
% INPUT
%   a lattice basis B = {b0, b1, ... , bn} in Z^m
%   a parameter £_ with 0.25 < £_ < 1, most commonly £_ = 0.75 (3/4)
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
return
