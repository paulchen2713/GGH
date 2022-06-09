%
% GGH public key cryptography with randomized v, w, r, with LLL reduction,
% and with no message input limits
%   ref. p.410 Example 7.36
%
clear;
clc;
%
fprintf('GGH public key cryptography with randomized v, w, r, with LLL reduction, and with no message input limits. \n\n');
%
% theoretically, n is the bigger the better, but here for our laptop's
% computation efficiency, we take n less or equal to 12, n <= 12
n = 10;
H_delta_v = 0.97^n;      % H(v) > 0.97^n
H_delta_w_u = 10 ^ (-4); % H(w) < 10^(-4), upper bound?
H_delta_w_d = 10 ^ (-5); % 10^(-5) < H)w), lower bound?
w_norm_delta = 10 ^(-4); % tolerance of the calculation error
power_v = 8;
power_w = 1;
power_r = 5;
% 
% random generation of v, and Hadamard ration H(v) has to be large enough
%
% iter = 0;
H_v = 0.00001;
while H_v < H_delta_v
    v = rand(n, n);
    v = round((10^power_v) * v) - 5 * 10^(power_v - 1);
    H_v = abs(det(v));
    for in = 1 : n
        H_v = H_v / norm(v(in, :));
    end
    H_v = H_v ^ (1/n);
%     iter = iter + 1;
%     fprintf('iter %d:  H_v = %f\n', iter, H_v);
%     display(H_v);
end
fprintf('H_v and v before LLL reduction\n');
fprintf('H_v = %f\n', H_v);
for in = 1 : n
    if in == 1
        fprintf('v = [ %d  %d  %d  %d  %d  %d  %d  %d  %d  %d;\n', v(in, :));
    elseif in > 1 && in < n
        fprintf('      %d  %d  %d  %d  %d  %d  %d  %d  %d  %d;\n', v(in, :));
    elseif in == n
        fprintf('      %d  %d  %d  %d  %d  %d  %d  %d  %d  %d ]\n\n', v(in, :));
    end
end
% display(v);
%
% LLL process
%
v = LLL(v);
H_v = abs(det(v));
for in = 1 : n
    H_v = H_v / norm(v(in, :));
end
H_v = H_v ^ (1/n);
fprintf('H_v and v after LLL reduction\n');
fprintf('H_v = %f\n', H_v);
for in = 1 : n
    if in == 1
        fprintf('v = [ %d  %d  %d  %d  %d  %d  %d  %d  %d  %d;\n', v(in, :));
    elseif in > 1 && in < n
        fprintf('      %d  %d  %d  %d  %d  %d  %d  %d  %d  %d;\n', v(in, :));
    elseif in == n
        fprintf('      %d  %d  %d  %d  %d  %d  %d  %d  %d  %d ]\n\n', v(in, :));
    end
end
% display(v);
% 
% random generation of u, and determint det(u) has to be +1 or -1
% 
% random generation of basis w = u * v, and H(w) has to be small enough, so
% that its inverse element can be calculate correctly
%
iter = 0;
w_norm = 1;
while w_norm > w_norm_delta
    H_w = H_delta_w_u;
    while H_w >= H_delta_w_u || H_w <= H_delta_w_d
        dd = rand(1, n);
        u = zeros(n, n);
        % randomly generate a diagonal matrix u
        % and det(u) has to be +1 or -1
        for in = 1 : n
            if dd(in) > 0.5
                u(in, in) = 1;
            else
                u(in, in) = -1;
            end
        end
        %
        % row operation
        %
        p1 = ceil(n * rand(n^2, 2));
        p2 = round((10 ^ power_w) * rand(1, n^2)) - 5 * 10^(power_w - 1);
        for in = 1 : n^2
            if p1(in, 1) ~= p1(in, 2)
                u(p1(in), :) = u(p1(in, 1), :) + p2(in) * u(p1(in, 2), :);
            else
                u(p1(in), :) = u(p1(in, 1), :) + p2(in) * u(mod(p1(in, 2) + 1, n) + 1, :);
            end
        end
        %
        % check whether u is an orthogonal matrix or not
        %
        det_u = det(u);
        fprintf('check whether u is an orthogonal matrix: ');
        fprintf(' det(u) = %f\n', det_u);
        %
        % compute a simple, but bad basis w = u * v, then
        % compute Hadamard ration H(w)
        %
        w = u * v;
        H_w = abs(det(w));
        for in = 1 : n
            H_w = H_w / norm(w(in, :));
        end
        H_w = H_w ^ (1/n);
        iter = iter + 1;
        fprintf('iter %d:  H_w = %f\n\n', iter, H_w);
%         display(H_w);
    end
    check_w = w * inv(w);
    w_norm = norm(check_w - eye(n));
    fprintf('After iter %d: w_norm = %f\n\n', iter, w_norm);
%     display(w_norm);
end
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
message1 = ' 1.Introduction to post-quantum cryptography';
message2 = ' 2.Quantum computing';
message3 = ' 3.Hash-based Digital Signature Schemes';
message4 = ' 4.Code-based cryptography';
message5 = ' 5.Lattice-based Cryptography';
message6 = ' 6.Multivariate Public Key Cryptography';
%
message_total = strcat(message1, message2, message3, message4, message5, message6);
message_len = length(message_total);
mN = ceil(message_len / n);
mR = mod(message_len, n);
for im = 1 : mN
    if im == mN
        message = message_total((im - 1)*n + 1 : message_len);
        if mR ~= 0
            for ir = mR+1 : n
                message(ir) = ' ';
            end
        end
    else
        message = message_total((im - 1)*n + 1 : im*n);
    end
    m = double(message);
    % 
    % random generation of r
    % 
    r = rand(1, n); % r = [ ... ]
    r = round((10 ^ power_r) * r) - 5 * 10^(power_r - 1);
    % 
    % compute the ciphertext e = m * w + r
    % 
    e(im, :) = m * w + r;
end
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Decryption
%
% use Babai's algorithm to compute the vector vv that belongs to thelattice to e
%
mN = size(e, 1);
message_r = zeros(1, 0);
for im = 1 : mN
    e_vv_coef = e(im, :) * inv(v);
    e_vv_coef = round(e_vv_coef);
    %
    vv = e_vv_coef * v;
    %
    % compute vv *w^(-1) to recover m
    % 
    m_r = vv * inv(w);
    m_r = round(m_r);
    %
    message_r = [message_r, m_r];
end
message_r = char(message_r);
%
%
fprintf('the original message is: %s\n', message_total);
fprintf('the recovery message is: %s\n\n', message_r);
%
%


