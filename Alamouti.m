% Mô phỏng STBC Alamouti 2x1 với MRC và ML
clear all; close all; clc;

frmLen = 1000; % chiều dài frame (số ký hiệu)
numPackets = 1000; % số lượng packet
EbN0 = 0:2:20; % giá trị Eb/N0 từ 0 đến 20 dB, bước 2 dB
N = 2; % số lượng anten Tx
M = 1; % số lượng anten thu Rx
P = 4; % độ dài điều chế QPSK

% Khởi tạo mảng lưu trữ BER
BERMRC = zeros(1, length(EbN0)); % BER cho MRC
BERML = zeros(1, length(EbN0)); % BER cho ML
errorMRC = zeros(1, numPackets); % Lỗi cho MRC
errorML = zeros(1, numPackets); % Lỗi cho ML

% Vòng lặp cho Eb/N0 điểm
for idx = 1:length(EbN0)
    for packetIdx = 1:numPackets
        % Tạo dữ liệu QPSK
        data = randi([0 P-1], frmLen, 2); % Dữ liệu QPSK (frmLen x 2)
        X = pskmod(data, P, 0); % Điều chế QPSK

        % Mã hóa Alamouti
        x1 = X(:,1); % Ký hiệu x1 (frmLen x 1)
        x2 = X(:,2); % Ký hiệu x2 (frmLen x 1)

        % Kênh Rayleigh fading
        Hr = (randn(frmLen, N, M) + 1i*randn(frmLen, N, M)) / sqrt(2);
        h1 = Hr(:,1,1); % Kênh từ anten 1
        h2 = Hr(:,2,1); % Kênh từ anten 2

        % Truyền tín hiệu
        s1 = h1.*x1 + h2.*x2; % Tín hiệu tại khe thời gian 1
        s2 = -h1.*conj(x2) + h2.*conj(x1); % Tín hiệu tại khe thời gian 2
        r1 = awgn(s1, EbN0(idx) + 3, 'measured'); % Thêm nhiễu (QPSK: Es/N0 = Eb/N0 + 3 dB)
        r2 = awgn(s2, EbN0(idx) + 3, 'measured');

        % Giải mã MRC
        norm = abs(h1).^2 + abs(h2).^2;
        x1Hat_MRC = (conj(h1).*r1 + h2.*conj(r2))./norm;
        x2Hat_MRC = (conj(h2).*r1 - h1.*conj(r2))./norm;

        % Gộp tín hiệu MRC
        decodedSymbols_MRC = [x1Hat_MRC x2Hat_MRC];
        decodedData_MRC = pskdemod(decodedSymbols_MRC, P, 0);

        % Tính lỗi cho MRC
        errorMRC(packetIdx) = biterr(decodedData_MRC, data);

    end

    % Tính BER (QPSK: 2 bit/ký hiệu, 2 cột dữ liệu)
    BERMRC(idx) = sum(errorMRC) / (numPackets * frmLen * 2 * 2);
    fprintf('Eb/N0 = %d dB, BER = %e\n', EbN0(idx), BERMRC(idx));
end

% Vẽ đồ thị
h = figure;
semilogy(EbN0, BERMRC, 'r*-', 'LineWidth', 1.5, 'DisplayName', 'Alamouti 2x1');
grid on;
set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-5 1]);
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('Hệ thống Alamouti 2x1 MRC');
legend('show');
set(h, 'color', [1 1 1]);
set(h, 'NumberTitle', 'off');
set(h, 'renderer', 'zbuffer');