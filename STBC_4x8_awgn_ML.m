% Mô phỏng STBC cho hệ thống MIMO 4x4 với giải mã Maximum Likelihood
clear all; close all; clc;

frmLen = 100; % chiều dài frame
numPackets = 1000; % số lượng packet
EbN0 = 0:2:10; % giá trị Eb/N0 từ 0 đến 20 dB, bước 2 dB
N = 4; % số lượng anten Tx
M = 4; % số lượng anten thu Rx
P = 4; % độ dài điều chế QPSK

% Khởi tạo
BER_ML = zeros(1, length(EbN0)); % Lưu BER cho ML
QPSK_const = pskmod(0:P-1, P, 0); % Tập hợp các điểm QPSK

% Vòng lặp cho EbN0 điểm
for idx = 1:length(EbN0)
    error_ML = zeros(1, numPackets); % Lưu lỗi cho mỗi packet
    
    % Vòng lặp cho số lượng packet
    for packetIdx = 1:numPackets
        data = randi([0 P-1], frmLen, 4); % Dữ liệu QPSK (4 cột)
        X = pskmod(data, P, 0); % Điều chế QPSK
        
        % Ma trận tín hiệu phát theo cấu trúc STBC
        x1 = X(:, 1:4);
        x2(:, 1) = -X(:, 2); x2(:, 2) = X(:, 1); x2(:, 3) = -X(:, 4); x2(:, 4) = X(:, 3);
        x3(:, 1) = -X(:, 3); x3(:, 2) = X(:, 4); x3(:, 3) = X(:, 1); x3(:, 4) = -X(:, 2);
        x4(:, 1) = -X(:, 4); x4(:, 2) = -X(:, 3); x4(:, 3) = X(:, 2); x4(:, 4) = X(:, 1);
        
        % Ma trận kênh
        Hr = (randn(frmLen, N, M) + 1i*randn(frmLen, N, M)) / sqrt(2);
        
        % Tín hiệu thu tại các anten
        r = zeros(frmLen, M, 4); % Lưu tín hiệu thu qua 4 thời điểm
        for n = 1:M
            %ma tran phat
            H = reshape(Hr(:,:,n),frmLen,N);
            %tin hieu tai anten thu
            r1(:,n)=awgn(sum(H.*x1,2)/sqrt(N),EbN0(idx));
            r2(:,n)=awgn(sum(H.*x2,2)/sqrt(N),EbN0(idx));
            r3(:,n)=awgn(sum(H.*x3,2)/sqrt(N),EbN0(idx));
            r4(:,n)=awgn(sum(H.*x4,2)/sqrt(N),EbN0(idx));

            z1(:,n) = r1(:,n).*conj(H(:,1)) + r2(:,n).*conj(H(:,2)) + r3(:,n).*conj(H(:,3)) + r4(:,n).*conj(H(:,4));

            z2(:,n) = r1(:,n).*conj(H(:,2)) - r2(:,n).*conj(H(:,1)) - r3(:,n).*conj(H(:,4)) + r4(:,n).*conj(H(:,3));

            z3(:,n) = r1(:,n).*conj(H(:,3)) + r2(:,n).*conj(H(:,4)) - r3(:,n).*conj(H(:,1)) - r4(:,n).*conj(H(:,2));

            z4(:,n) = r1(:,n).*conj(H(:,4)) - r2(:,n).*conj(H(:,3)) + r3(:,n).*conj(H(:,2)) - r4(:,n).*conj(H(:,1));
        end
        r = [sum(z1,2) sum(z2,2) sum(z3,2) sum(z4,2)];

        demod22 = pskdemod(r, P, 0);
        error22(packetIdx) = biterr(demod22,data);
    end
        BER22(idx) = sum(error22)/(numPackets*frmLen);
        fprintf('Eb/N0 = %d dB, BER = = %e\n', EbN0(idx), BER22(idx));
        semilogy(EbN0(1:idx),BER22(1:idx),'r*-');
end
% Vẽ đồ thị
h = figure;
semilogy(EbN0, BER22, 'r*-', 'LineWidth', 1.5, 'DisplayName', 'MIMO 4x4');
grid on;
set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-5 1]);
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('Hệ thống MIM0 4x4');
legend('show');
set(h, 'color', [1 1 1]);
set(h, 'NumberTitle', 'off');
set(h, 'renderer', 'zbuffer');