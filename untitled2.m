clear all; close all; clc;
frmLen = 2; % Chiều dài frame
N = 2; % 2 anten phát
M = 2; % 2 anten thu
snr = 10; % SNR = 10 dB

% Tạo tín hiệu truyền (giả sử đã mã hóa STBC, dùng Alamouti cho đơn giản)
tx2 = [1+0i  0+1i; ...
       0-1i  1+0i]; % [frmLen x N] = [2 x 2]

% Tạo ma trận kênh H
H = (randn(frmLen, N, M) + 1j*randn(frmLen, N, M))/sqrt(2);

% Tính tín hiệu nhận tại anten thu thứ m
rxSymbols = zeros(frmLen, M);
for m = 1:M
    % Bước 1: Tín hiệu qua kênh
    signal_after_channel = H(:, :, m).*tx2;
    disp("Tín hiệu qua kênh ")
    disp(signal_after_channel)
    % Bước 2: Tổng tín hiệu từ các anten phát
    signal_sum = sum(signal_after_channel, 2);
    % Bước 3: Chuẩn hóa công suất
    signal_normalized = signal_sum/sqrt(N);
    % Bước 4: Thêm nhiễu
    rxSymbols(:, m) = awgn(signal_normalized, snr);
end

% Hiển thị kết quả
disp('Tín hiệu truyền (tx2):');
disp(tx2);
disp('Ma trận kênh H:');
disp(H);
disp('Tín hiệu nhận (rxSymbols):');
disp(rxSymbols);