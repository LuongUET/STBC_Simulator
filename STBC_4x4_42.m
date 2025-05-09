clear all; close all; clc;
frmLen = 1e4; % chiều dài frame
numPackets = 100; % số lượng packet
EbN0 = 0:0.5:10; % giá trị Eb/N0 từ 0 đến 20 dB, bước 2 dB
N = 4; % số lượng anten Tx
M = 4; % số lượng anten thu Rx
P = 4; % độ dài điều chế QPSK

% khởi tạo
tx2 = zeros(frmLen,N); H = zeros(frmLen,N,M);
r21 = zeros(frmLen,2); r12 = zeros(frmLen,2);r22 = zeros(frmLen,2,2);
z21 = zeros(frmLen,1); z21_1 = zeros(frmLen/N,1); z21_2 = z21_1;
z12 = zeros(frmLen,M); z22_1 = zeros(frmLen/N,1); z22_2 = z22_1; 
z22 = zeros(frmLen,2); zr22 = zeros(frmLen,1);
error11 = zeros(1,numPackets); BER11 = zeros(1,length(EbN0));
error21 = error11; BER21 = BER11; error12 = error11; BER12 = BER11;
BERthy2 = BER11;

% vòng lặp cho EbN0 điểm
for idx = 1:length(EbN0)
    % vòng lặp cho số lượng packet
    for packetIdx = 1:numPackets
        data = randi([0 P-1], frmLen, 4); % Dữ liệu QPSK (4 cột)
        X = pskmod(data, P,0); % điều chế QPSK
        % khởi tạo figure cho kết quả BER
        x0 = X(:,1);
        x1 = X(:,1:4);
        x2(:,1) = -X(:,2); x2(:,2) = X(:,1); x2(:,3) = -X(:,4); x2(:,4) = X(:,3);
        x3(:,1) = -X(:,3); x3(:,2) = X(:,4); x3(:,3) = X(:,1); x3(:,4) = -X(:,2);
        x4(:,1) = -X(:,4); x4(:,2) = -X(:,3); x4(:,3) = X(:,2); x4(:,4)= X(:,1);
   
        x5 = conj(x1);
        x6 = conj(x2);
        x7 = conj(x3);
        x8 = conj(x4);

        Hr = (randn(frmLen, N, 4, M) + 1i*randn(frmLen, N, 4, M))/sqrt(2);
        z1 = zeros(frmLen,4,M);
        z2 = zeros(frmLen,4,M);
        z3 = zeros(frmLen,4,M);
        z4 = zeros(frmLen,4,M);
        for mRx = 1:M
            for n=1:4
                %ma tran phat
                 H = reshape(Hr(:,:,n,mRx), frmLen, N);
                %tin hieu tai anten thu
                r1(:,n)=awgn(sum(H.*x1,2)/sqrt(N),EbN0(idx));
                r2(:,n)=awgn(sum(H.*x2,2)/sqrt(N),EbN0(idx));
                r3(:,n)=awgn(sum(H.*x3,2)/sqrt(N),EbN0(idx));
                r4(:,n)=awgn(sum(H.*x4,2)/sqrt(N),EbN0(idx));
                r5(:,n)=awgn(sum(H.*x5,2)/sqrt(N),EbN0(idx));
                r6(:,n)=awgn(sum(H.*x6,2)/sqrt(N),EbN0(idx));
                r7(:,n)=awgn(sum(H.*x7,2)/sqrt(N),EbN0(idx));
                r8(:,n)=awgn(sum(H.*x8,2)/sqrt(N),EbN0(idx));
                
                z1_1(:,n) = r1(:,n).*conj(H(:,1)) + r2(:,n).*conj(H(:,2)) + r3(:,n).*conj(H(:,3)) + r4(:,n).*conj(H(:,4));
                z1_2(:,n) = conj(r5(:,n)).*H(:,1) + conj(r6(:,n)).*H(:,2) + conj(r7(:,n)).*H(:,3) + conj(r8(:,n)).*H(:,4);
                z1(:,n,mRx) = z1_1(:,n) + z1_2(:,n);
    
                z2_1(:,n) = r1(:,n).*conj(H(:,2)) - r2(:,n).*conj(H(:,1)) - r3(:,n).*conj(H(:,4)) + r4(:,n).*conj(H(:,3));
                z2_2(:,n) = conj(r5(:,n)).*H(:,2) - conj(r6(:,n)).*H(:,1) - conj(r7(:,n)).*H(:,4) + conj(r8(:,n)).*H(:,3);
                 z2(:,n,mRx) = z2_1(:,n) + z2_2(:,n);
    
                z3_1(:,n) = r1(:,n).*conj(H(:,3)) + r2(:,n).*conj(H(:,4)) - r3(:,n).*conj(H(:,1)) - r4(:,n).*conj(H(:,2));
                z3_2(:,n) = conj(r5(:,n)).*H(:,3) + conj(r6(:,n)).*H(:,4) - conj(r7(:,n)).*H(:,1) - conj(r8(:,n)).*H(:,2);
                z3(:,n,mRx) = z3_1(:,n) + z3_2(:,n);
    
                z4_1(:,n) = r1(:,n).*conj(H(:,4)) - r2(:,n).*conj(H(:,3)) + r3(:,n).*conj(H(:,2)) - r4(:,n).*conj(H(:,1));
                z4_2(:,n) = conj(r5(:,n)).*H(:,4) - conj(r6(:,n)).*H(:,3) + conj(r7(:,n)).*H(:,2) - conj(r8(:,n)).*H(:,1);
                z4(:,n,mRx) = z4_1(:,n) + z4_2(:,n);    
                
            end
        end
        r = [sum(sum(z1,3),2)  sum(sum(z2,3),2)  sum(sum(z3,3),2)  sum(sum(z4,3),2)];  
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
set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-1 1e-5]);
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('Hệ thống MIM0 4x4');
legend('show');
set(h, 'color', [1 1 1]);
set(h, 'NumberTitle', 'off');
set(h, 'renderer', 'zbuffer');