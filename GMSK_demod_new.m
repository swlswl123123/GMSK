function [decode_out] = GMSK_demod_new(signal_recv, decode, i, g)
    %myFun - Description
    %
    % Syntax: output = myFun(input)
    %
    % Long description
    dif_min = 100;
    decode_out = decode;
    num_bits_pulse = 304;
    path = zeros(2^3, 3);
    f_c = 240e6;
    
    path(1,:) = [0,0,0];
    path(2,:) = [0,0,1];
    path(3,:) = [0,1,0];
    path(4,:) = [0,1,1];
    path(5,:) = [1,0,0];
    path(6,:) = [1,0,1];
    path(7,:) = [1,1,0];
    path(8,:) = [1,1,1];
    
    if i == 1
        path = [repmat(zeros(1,2),8,1),path];
    elseif i == 2
        path = [repmat(zeros(1,1),8,1),repmat(decode(i-1),8,1),path];
    else
        path = [repmat(decode(i-2:i-1),8,1), path];
    end
    
    path = path * 2 - 1;
    
    for j = 1:8
        phi_last_tmp = 0;
        bit_5 = path(j,:);
        [~, I_sig, Q_sig, ~] = GMSK(bit_5, f_c, phi_last_tmp, g);
        g_path(j,:) = complex(I_sig, Q_sig);
    
        g_path_sample = g_path(j, 1:16:end);
        signal_recv_ang = angle(signal_recv);
        g_path_ang = angle(g_path_sample);
    
        % signal_recv_ang = signal_recv_ang + (g_path_dif_ang(1) - signal_recv_ang(1));
        % signal_recv_ang(signal_recv_ang >= 2*pi) = signal_recv_ang(signal_recv_ang >= 2*pi) - 2*pi;
        % signal_recv_ang(signal_recv_ang <= -2*pi) = signal_recv_ang(signal_recv_ang <= -2*pi) + 2*pi;
    
        if i == 5
            figure
            plot(g_path_ang);
            hold on;
            plot(signal_recv_ang);
            hold off;
        end
    
        g_path_ang_dif = zeros(1,3);
        signal_recv_ang_dif = zeros(1,3);
    
        for m = 1:3
            g_path_ang_dif(m) = g_path_ang(m) - g_path_ang(m+1);
            if  g_path_ang_dif(m) > pi
                g_path_ang_dif(m) =  g_path_ang_dif(m) - 2*pi;
            elseif  g_path_ang_dif(m) < -pi
                g_path_dif_dif(m) =  g_path_ang_dif(m) + 2*pi;
            end
            signal_recv_ang_dif(m) = signal_recv_ang(m) - signal_recv_ang(m+1);
            if signal_recv_ang_dif(m) > pi
                signal_recv_ang_dif(m) = signal_recv_ang_dif(m) - 2*pi;
            elseif signal_recv_ang_dif(m) < -pi
                signal_recv_ang_dif(m) = signal_recv_ang_dif(m) + 2*pi;
            end
        end

        % if i == 5
        %     signal_recv_ang_dif
        %     dif(j) = (sum(signal_recv_ang_dif) - sum(g_path_ang_dif)).^2
        % end
    
        % dif(j) = sum((g_path_dif_ang - signal_recv_ang).^2)
        dif(j) = (sum(signal_recv_ang_dif) - sum(g_path_ang_dif)).^2;
        if dif(j) < dif_min
            dif_min = dif(j);
            decode_out(i) = (path(j,3) + 1)/2;
        end
    end
    
    
    
    % for n = 1:16
    %     if path(n, 3) == -1
    %         dif_zero = dif_zero + dif(n);
    %     else
    %         dif_one = dif_one + dif(n);
    %     end
    % end
    
    % % dif_one
    % % dif_zero
    
    % if dif_one < dif_zero
    %     decode_out(i) = 1;
    % else
    %     decode_out(i) = 0;
    % end
        
    end