% Mô phỏng STBC cho hệ thống MIMO 4x4
clear all; close all; clc;

% Cấu hình hệ thống
frmLen = 100;        % chiều dài frame
numPackets = 1000;   % số lượng packets
EbN0 = 0:2:20;       % giá trị Eb/N0 từ 0 đến 20 dB
N = 4;               % số lượng anten phát Tx
M = 4;               % số lượng anten thu Rx
P = 4;               % QPSK (4-QAM)

% Khởi tạo biến để lưu kết quả
ber = zeros(1, length(EbN0));

% Vòng lặp qua các giá trị Eb/N0
for idx = 1:length(EbN0)
    snr = EbN0(idx);
    fprintf('SNR = %d dB\n', snr);
    
    % Khởi tạo biến đếm lỗi
    totalErrors = 0;
    totalBits = 0;
    
    % Vòng lặp qua các packet
    for packetIdx = 1:numPackets
        % Tạo dữ liệu ngẫu nhiên
        data = randi([0 P-1], frmLen, 4); % Dữ liệu QPSK (4 cột)

        % Điều chế QPSK
        txSymbols = pskmod(data, P, 0); % Bỏ 'SymbolOrder', 'gray' vì không được hỗ trợ
        
        % Mã hóa STBC với ma trận X_4
        x1 = txSymbols(:, 1); x2 = txSymbols(:, 2); x3 = txSymbols(:, 3); x4 = txSymbols(:, 4);
        tx2 = zeros(frmLen, N);
        for t = 1:4:frmLen
            if t+3 <= frmLen
                % Tạo ma trận X_4 cho 4 ký hiệu tại thời điểm t
                 X4 = [x1(t) -x2(t) -x3(t) -x4(t); ...
                       x2(t)  x1(t)  x4(t) -x3(t); ...
                       x3(t) -x4(t)  x1(t)  x2(t); ...
                       x4(t)  x3(t) -x2(t)  x1(t)];
                    % Gán vào tx2
                 tx2(t:t+3, :) = X4;
            end
        end
        
        % Kênh truyền (Rayleigh fading)
        H = (randn(frmLen, N, M) + 1j*randn(frmLen, N, M))/sqrt(2);
        
        % Tín hiệu thu
        rxSymbols = zeros(frmLen, M);
        for m = 1:M
            rxSymbols(:, m) = awgn(sum(H(:, :, m).*tx2, 2)/sqrt(N), snr);
        end
        
        % Giải mã STBC
        z1 = zeros(frmLen, 1); z2 = z1; z3 = z1; z4 = z1;
        for t = 1:4:frmLen
            if t+3 <= frmLen
                y = rxSymbols(t:t+3, :);
                H_t = H(t:t+3, :, :);
                z1(t:t+3) = sum(squeeze(conj(H_t(:, 1, :))).*y, 2);
                z2(t:t+3) = sum(squeeze(conj(H_t(:, 2, :))).*y, 2);
                z3(t:t+3) = sum(squeeze(conj(H_t(:, 3, :))).*y, 2);
                z4(t:t+3) = sum(squeeze(conj(H_t(:, 4, :))).*y, 2);
            end
        end
        
        % Giải điều chế
        rxSymbolsDecoded = [z1 z2 z3 z4];
        rxBits = pskdemod(rxSymbolsDecoded, P, 0); % Bỏ 'SymbolOrder', 'gray'
        
        % Tính lỗi bit
        errors = sum(sum(abs(rxBits - data)));
        totalErrors = totalErrors + errors;
        totalBits = totalBits + numel(data);
    end
    
    % Tính BER
    ber(idx) = totalErrors / totalBits;
end

% Vẽ kết quả
figure;
semilogy(EbN0, ber, 'r*-', 'LineWidth', 2, 'DisplayName', 'BER (STBC 4x4)');
grid on;
xlabel('Eb/N0 (dB)');
ylabel('BER');
title('BER vs Eb/N0 for MIMO 4x4 with STBC');
legend('show');