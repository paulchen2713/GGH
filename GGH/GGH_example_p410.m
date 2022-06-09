%
% GGH public key cryptography
%   ref. p.410 Example 7.36
%
clear;
clc;
%
fprintf('GGH public key cryptography Example 7.36, p.410 \n\n');
%
% theoretically, n is the bigger the better, but here for our laptop's
% computation efficiency, we take n less or equal to 12, n <= 12
%
n = 3;
v = [-97 19 19; -36 30 86; -184 -64 78];
%
% compute the Hadamard ratio of v
%
H_v = abs(det(v));
for in = 1 : n
    H_v = H_v / norm(v(in, :));
end
H_v = H_v ^ (1/3);
fprintf('the Hadamard ratio H_v = %f\n\n', H_v);
%
%
u = [4327 -15447 23454; 3297 -11770 17871; 5464 -19506 29617];
%
% check whether u is an orthogonal matrix
%
fprintf('check whether u is an orthogonal matrix: ');
fprintf(' det(u) = %f\n\n', det(u));
% display(det(u));
%
% compute a bad basis w = u * v, then 
% compute the Hadamard ratio of w
%
w = u * v;
fprintf('the basis w = u * v\n');
for in = 1 : n
    if in == 1
        fprintf('w = [ %d  %d  %d;\n', w(in, :));
    elseif in > 1 && in < n
        fprintf('      %d  %d  %d;\n', w(in, :));
    elseif in == n
        fprintf('      %d  %d  %d ]\n\n', w(in, :));
    end
end
% display(w);
H_w = abs(det(w));
for in = 1 : n
    H_w = H_w / norm(w(in, :));
end
H_w = H_w ^ (1/3);
fprintf('the Hadamard ratio H_w = %f\n\n', H_w);
% display(H_w);
%
% private key: v, u
% public key:  w
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Encryption
%
% message: m
%
m = [86 -35 -32];
%
r = [-4 -3 2];
%
%  compute the ciphertext e = m * w + r
%
e = m * w + r;
fprintf('the ciphertext e = m * w + r\n');
fprintf('e = [%d %d %d]\n\n', e);
% display(e);
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Decryption
%
% use Babai's algorithm to compute the vector vv that belongs to thelattice to e
% 
e_vv_coef = e * inv(v);
% display(e_vv_coef);
e_vv_coef = round(e_vv_coef);
fprintf('the coefficient e_vv_coef = e * inv(v)\n');
fprintf('e_vv_coef = [%d %d %d]\n\n', e_vv_coef);
% display(e_vv_coef);
%
vv = e_vv_coef * v;
fprintf('the vector vv = e_vv_coef * v\n');
fprintf('vv = [%d %d %d]\n\n', vv);
% display(vv);
%
% compute vv *w^(-1)
% 
m_r = vv * inv(w);
fprintf('m_r = [%f %f %f]\n\n', m_r);
% display(m_r);
%
%
fprintf('the original message m   is: [%d %d %d]\n', m);
% display(m);
fprintf('the recovery message m_r is: [%f %f %f]\n\n', m_r);
% display(m_r);
%
%

