clear all; 
close all;

bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间

f_IF = 240e6; %射频频率
fs_IF = 1024e6;  % 射频、中频频信号采样速率
fs_BB = 128e6;  % 基带信号采样速率
num_bits_pulse = 304; 
oversamp_BB = T * fs_BB;  % 基带信号过采样速率
oversamp_IF = T * fs_IF;  % 射频、中频信号过采样速率
T_s_BB = 1/fs_BB;  % 基带采样间隔

load('lib/g_1024.mat');  % GMSK调制 g函数 

% generate code

I = randi(2,1,num_bits_pulse);
I = I - 1;

% coding
I = 2*I - 1;
bit_5 = zeros(1,5);
phi_last = 0;

for i = 1:num_bits_pulse
    if i == 1
        bit_5 = [0,0,I(i:i+2)];
    elseif i == 2
        bit_5 = [0,I(i:i+3)];
    elseif i == num_bits_pulse-1
        bit_5 = [I(i-2:i+1),0];
    elseif i == num_bits_pulse
        bit_5 = [I(i-2:i),0,0];
    else
        bit_5 = I(i-2:i+2);
    end

    [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, f_IF, phi_last, g);
    signal_trans_temp_BB((i-1)*oversamp_IF+1:(i)*oversamp_IF) = complex(I_sig, Q_sig);
    phi_all((i-1)*oversamp_IF+1:(i)*oversamp_IF) = phi_int;
    plot(mod(phi_all,2*pi))
end