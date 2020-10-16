function [decode_out, phi_last_out] = GMSK_demod(signal_recv_dif, decode, i, phi_last, g)
%myFun - Description
%
% Syntax: output = myFun(input)
%
% Long description
decode_out = decode;
path = zeros(2^4, 4);
g_path = zeros(16, 64*2);
g_path_dif = zeros(16, 64);
tendency = zeros(1, 16);
tendency_max = 0;
f_c = 240e6;

path(1,:) = [0,0,0,0];
path(2,:) = [0,0,0,1];
path(3,:) = [0,0,1,0];
path(4,:) = [0,0,1,1];
path(5,:) = [0,1,0,0];
path(6,:) = [0,1,0,1];
path(7,:) = [0,1,1,0];
path(8,:) = [0,1,1,1];
path(9,:) = [1,0,0,0];
path(10,:) = [1,0,0,1];
path(11,:) = [1,0,1,0];
path(12,:) = [1,0,1,1];
path(13,:) = [1,1,0,0];
path(14,:) = [1,1,0,1];
path(15,:) = [1,1,1,0];
path(16,:) = [1,1,1,1];

if i == 1
    path = [repmat(zeros(1,2),16,1),path];
elseif i == 2
    path = [repmat(zeros(1,1),16,1), repmat(decode(i-1),16,1), path];
else
    path = [repmat(decode(i-2:i-1),16,1), path];
end

path = path * 2 - 1;

for j = 1:16
    phi_last_tmp = phi_last;
    phi_all = zeros(1,64*2);
    for k = 1:2
        bit_5 = path(j,k:k+4);
        [phi_last_tmp, I_sig, Q_sig, phi_all(1+(k-1)*64:k*64)] = GMSK(bit_5, f_c, phi_last_tmp, g);
        g_path(j, 1+(k-1)*64 : k*64) = complex(I_sig, Q_sig);
    end
    figure
    % plot(phi_all)
    g_path_dif(j,:) = g_path(j, 1:64).*conj(g_path(j, 65:128));
    plot(angle(g_path_dif(j,:)))
    hold on;
    plot(angle(signal_recv_dif));
    hold off;

    g_path_dif_tendency = zeros(1,63);
    signal_recv_dif_tendency = zeros(1,63);

    for m = 1:63
        g_path_dif_tendency(m) = angle(g_path_dif(j,m)) - angle(g_path_dif(j,m+1));
        signal_recv_dif_tendency(m) = angle(signal_recv_dif(m)) - angle(signal_recv_dif(m+1));
    end

    tendency(j) = sum(g_path_dif_tendency.*signal_recv_dif_tendency)
    if tendency(j) > tendency_max
        tendency_max = tendency(j);
        phi_last_out = phi_last_tmp;
        decode_out(i) = (path(j,3) + 1)/2;
    end
end



% for n = 1:16
%     if path(n, 3) == -1
%         tendency_zero = tendency_zero + tendency(n);
%     else
%         tendency_one = tendency_one + tendency(n);
%     end
% end

% % tendency_one
% % tendency_zero

% if tendency_one < tendency_zero
%     decode_out(i) = 1;
% else
%     decode_out(i) = 0;
% end
    
end