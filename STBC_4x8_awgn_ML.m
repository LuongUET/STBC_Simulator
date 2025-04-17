% Mô phỏng STBC cho hệ thống MIMO 4x4 với giải mã Maximum Likelihood
clear all; close all; clc;

frmLen = 100; % chiều dài frame
numPackets = 1; % số lượng packet
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
            H = reshape(Hr(:, :, n), frmLen, N);
            r(:, n, 1) = awgn(sum(H .* x1, 2) / sqrt(N), EbN0(idx));
            r(:, n, 2) = awgn(sum(H .* x2, 2) / sqrt(N), EbN0(idx));
            r(:, n, 3) = awgn(sum(H .* x3, 2) / sqrt(N), EbN0(idx));
            r(:, n, 4) = awgn(sum(H .* x4, 2) / sqrt(N), EbN0(idx));
        end
        
        % Giải mã Maximum Likelihood
        demod_ML = zeros(frmLen, 4); % Lưu dữ liệu giải mã
        for t = 1:frmLen
            % Tín hiệu thu tại thời điểm t
            r_t = squeeze(r(t, :, :)).'; % [4 x M], 4 thời điểm, M anten thu
            
            % Ma trận kênh tại thời điểm t
            H_t = squeeze(Hr(t, :, :)); % [N x M]
            
            min_dist = inf; % Khoảng cách nhỏ nhất
            best_sym = zeros(1, 4); % Ký hiệu tốt nhất
            
            % Duyệt qua tất cả tổ hợp của 4 ký hiệu QPSK
            for s1 = 1:P
                for s2 = 1:P
                    for s3 = 1:P
                        for s4 = 1:P
                            % Tín hiệu thử
                            x_try = [QPSK_const(s1), QPSK_const(s2), QPSK_const(s3), QPSK_const(s4)];
                            
                            % Tạo khối STBC cho tín hiệu thử
                            x1_try = x_try;
                            x2_try = [-x_try(2), x_try(1), -x_try(4), x_try(3)];
                            x3_try = [-x_try(3), x_try(4), x_try(1), -x_try(2)];
                            x4_try = [-x_try(4), -x_try(3), x_try(2), x_try(1)];
                            
                            % Tín hiệu dự đoán
                            r_pred = zeros(4, M);
                            for m = 1:M
                                r_pred(1, m) = sum(H_t(:, m) .* x1_try.') / sqrt(N);
                                r_pred(2, m) = sum(H_t(:, m) .* x2_try.') / sqrt(N);
                                r_pred(3, m) = sum(H_t(:, m) .* x3_try.') / sqrt(N);
                                r_pred(4, m) = sum(H_t(:, m) .* x4_try.') / sqrt(N);
                            end
                            
                            % Tính khoảng cách Euclidean
                            dist = norm(r_pred - r_t, 'fro')^2;
                            
                            % Cập nhật nếu khoảng cách nhỏ hơn
                            if dist < min_dist
                                min_dist = dist;
                                best_sym = x_try;
                            end
                        end
                    end
                end
            end
            
            % Giải điều chế ký hiệu tốt nhất
            demod_ML(t, :) = pskdemod(best_sym, P, 0);
        end
        
        % Tính lỗi bit
        error_ML(packetIdx) = biterr(demod_ML, data);
    end
    
    % Tính BER
    BER_ML(idx) = sum(error_ML) / (numPackets * frmLen * log2(P));
    
    % Vẽ đồ thị
    figure(1); grid on; hold on;
    semilogy(EbN0(1:idx), BER_ML(1:idx), 'b*-', 'DisplayName', 'ML Decoding');
    xlabel('Eb/N0 (dB)'); ylabel('BER');
    title('Hệ thống MIMO 4x4 với STBC');
    set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-5 1]);
    legend('show');
end