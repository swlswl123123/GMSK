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
load('lpf_coe.mat');  % 21个频点

% generate code

I_single = randi(2,1,num_bits_pulse);
I_single = I_single - 1;
I = 2*I_single - 1;
% I = [-1,-1,1,-1,-1,-1,1,1,1,1];

% coding
bit_5 = zeros(1,5);
phi_last = 0;

for i = 1:num_bits_pulse
    if i == 1
        bit_5 = [-1,-1,I(i:i+2)];
    elseif i == 2
        bit_5 = [-1,I(i-1:i+2)];
    elseif i == num_bits_pulse-1
        bit_5 = [I(i-2:i+1),-1];
    elseif i == num_bits_pulse
        bit_5 = [I(i-2:i),-1,-1];
    else
        bit_5 = I(i-2:i+2);
    end

    [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, f_IF, phi_last, g);
    signal_trans_BB((i-1)*oversamp_IF+1:(i)*oversamp_IF) = complex(I_sig, Q_sig);
    phi_all((i-1)*oversamp_IF+1:(i)*oversamp_IF) = phi_int;
end
% figure;
% plot(mod(phi_all,2*pi))

% 射频信号

t = linspace(0, num_bits_pulse*T, oversamp_IF*num_bits_pulse);
signal_trans_IF = signal_trans_BB .* exp(1i*2*pi*f_IF*t);
% signal_trans_IF = exp(1i*2*pi*f_IF*t);
% plot(t, signal_trans_IF)

% 接收端
% 加噪声
SNRdB = 8;
signal_recv_IF_noise = awgn(signal_trans_IF, SNRdB, 'measured');
% signal_recv_IF_noise = signal_trans_IF;
signal_recv_noise_IF_FFT = abs(fft(signal_recv_IF_noise));

fre = t./(num_bits_pulse*T)*fs_IF;
% plot(fre, signal_recv_noise_IF_FFT)
% hold on;

% 理想带通滤波

BandPass = zeros(size(fre));
BandPass(round((f_IF-70e6)/fs_IF*oversamp_IF*num_bits_pulse):round((f_IF+70e6)/fs_IF*oversamp_IF*num_bits_pulse)) = 1;
signal_recv_IF = ifft(BandPass.*fft(signal_recv_IF_noise));
% signal_recv_IF = signal_recv_IF_noise;

% plot(fre, abs(fft(signal_recv_IF)))

% 译码
decode = zeros(size(I));
% for i = 1:2:num_bits_pulse
for i = 1:2:3
    signal_recv_dif = signal_recv_IF(1+(i-1)*oversamp_IF:i*oversamp_IF).*conj(signal_recv_IF(1+i*oversamp_IF:(i+1)*oversamp_IF));
    decode = GMSK_demod(signal_recv_dif, decode, i, g);
end

error = I_single - decode;

error(error~=0) = 1;

error_rate = sum(error)/304




