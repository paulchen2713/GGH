% 
% Gram-Schmidt orthogonalized process function
% 
function BB = GramSchmidt(B)
B_size = size(B, 1);
BB = B;
for in = 2 : B_size
    for j = 1 : in - 1
        BB(in, :) = BB(in, :) - (BB(in, :) * BB(j, :)' / (BB(j, :) * BB(j, :)')) * BB(j, :);
    end
end
return
