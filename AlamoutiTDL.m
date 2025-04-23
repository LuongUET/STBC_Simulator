% Mô phỏng STBC Alamouti 2x1 với kênh TDL (5G)
clear all; close all; clc;

frmLen = 1000; % chiều dài frame (số ký hiệu)
numPackets = 1000; % số lượng packet
EbN0 = 0:2:20; % giá trị Eb/N0 từ 0 đến 20 dB, bước 2 dB
N = 2; % số anten Tx
M = 1; % số anten thu Rx
P = 4; % độ dài điều chế QPSK

% Tham số kênh TDL-A (theo 3GPP TR 38.901)
numTaps = 23; % Số tap của TDL-A
delays = [0, 30, 70, 90, 110, 190, 410, 50, 120, 200, 230, 500, 1600, ...
          2300, 2500, 1400, 2100, 2700, 2900, 3100, 3200, 3400, 3700]; % Độ trễ (ns)
power_dB = [0, -1, -2, -3, -8, -7.8, -7.5, -15.9, -6.2, -16.6, -13.4, -18.6, ...
            -20.8, -16, -19.2, -22.8, -17.6, -22.7, -24, -25.7, -23.9, -24.2, -25.2]; % Công suất (dB)
power = 10.^(power_dB/10); % Chuyển sang thang tuyến tính
power = power / sum(power); % Chuẩn hóa tổng công suất = 1

% Tham số Doppler
v = 3; % Tốc độ di chuyển (km/h)
v_ms = v * 1000 / 3600; % Chuyển sang m/s
fc = 3.5e9; % Tần số sóng mang (3.5 GHz)
c = 3e8; % Tốc độ ánh sáng
fd_max = v_ms * fc / c; % Tần số Doppler tối đa
Ts = 1e-6; % Chu kỳ lấy mẫu (1 μs, giả định)
dopplerRate = fd_max * Ts; % Tỷ lệ Doppler chuẩn hóa

% Khởi tạo mảng lưu trữ BER
BERMRC = zeros(1, length(EbN0)); % BER cho MRC
errorMRC = zeros(1, numPackets); % Lỗi cho MRC

% Vòng lặp cho Eb/N0
for idx = 1:length(EbN0)
    for packetIdx = 1:numPackets
        % Tạo dữ liệu QPSK
        data = randi([0 P-1], frmLen, 2); % Dữ liệu QPSK (frmLen x 2)
        X = pskmod(data, P, 0); % Điều chế QPSK

        % Mã hóa Alamouti
        x1 = X(:,1); % Ký hiệu x1 (frmLen x 1)
        x2 = X(:,2); % Ký hiệu x2 (frmLen x 1)

        % --- Tạo kênh TDL ---
        h1 = zeros(frmLen, 1); % Hệ số kênh tổng hợp cho anten 1
        h2 = zeros(frmLen, 1); % Hệ số kênh tổng hợp cho anten 2
        for l = 1:numTaps
            % Tạo fading Rayleigh cho mỗi tap
            hl = (randn(frmLen, 1) + 1i*randn(frmLen, 1)) / sqrt(2);
            % Áp dụng Doppler (mô phỏng đơn giản bằng lọc)
            hl = filter(1, [1 -0.99*exp(-2*pi*1i*dopplerRate)], hl);
            % Chuẩn hóa công suất tap
            hl = hl * sqrt(power(l));
            % Tổng hợp kênh cho anten 1
            h1 = h1 + hl;
            % Lặp lại cho anten 2 (kênh độc lập)
            hl = (randn(frmLen, 1) + 1i*randn(frmLen, 1)) / sqrt(2);
            hl = filter(1, [1 -0.99*exp(-2*pi*1i*dopplerRate)], hl);
            hl = hl * sqrt(power(l));
            h2 = h2 + hl;
        end
        % Chuẩn hóa công suất kênh tổng hợp
        h1 = h1 / sqrt(mean(abs(h1).^2)) * sqrt(sum(power));
        h2 = h2 / sqrt(mean(abs(h2).^2)) * sqrt(sum(power));

        % Truyền tín hiệu
        s1 = h1.*x1 + h2.*x2; % Tín hiệu tại khe thời gian 1
        s2 = -h1.*conj(x2) + h2.*conj(x1); % Tín hiệu tại khe thời gian 2
        r1 = awgn(s1, EbN0(idx) + 3, 'measured'); % Thêm nhiễu
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
semilogy(EbN0, BERMRC, 'r*-', 'LineWidth', 1.5, 'DisplayName', 'Alamouti 2x1 TDL-A');
grid on;
set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-5 1]);
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('Hệ thống Alamouti 2x1 với kênh TDL-A');
legend('show');
set(h, 'color', [1 1 1]);
set(h, 'NumberTitle', 'off');
set(h, 'renderer', 'zbuffer');