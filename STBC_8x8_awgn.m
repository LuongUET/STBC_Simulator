clc; clear;

N = 8; % Số anten phát
M = 8; % Số anten thu
P = 4; % QPSK
EbN0 = 0:2:20;
frmLen = 8; % Số frame để tính BER
numPackets = 1000; % số lượng packet
BER = zeros(1, length(EbN0));

for idx = 1:length(EbN0)
    errors = 0; total = 0;
    snr = EbN0(idx);
    
    for packetIdx = 1:numPackets
        %Tạo dữ liệu và điều chế QPSK
        data = randi([0 P-1], frmLen, 8); % 8 symbol
        x= pskmod(data, P, 0);  % QPSK
        
        % Time slot 1–8
        X_coded(1,:)  = [x(1) x(2) x(3) x(4) -conj(x(5)) -conj(x(6)) -conj(x(7)) -conj(x(8))];
        X_coded(2,:)  = [x(4) x(1) x(2) x(3) -conj(x(8)) -conj(x(5)) -conj(x(6)) -conj(x(7))];
        X_coded(3,:)  = [x(3) x(4) x(1) x(2) -conj(x(7)) -conj(x(8)) -conj(x(5)) -conj(x(6))];
        X_coded(4,:)  = [x(2) x(3) x(4) x(1) -conj(x(6)) -conj(x(7)) -conj(x(8)) -conj(x(5))];
        X_coded(5,:)  = [x(5) x(6) x(7) x(8) -conj(x(1)) -conj(x(2)) -conj(x(3)) -conj(x(4))];
        X_coded(6,:)  = [x(8) x(5) x(6) x(7) -conj(x(4)) -conj(x(1)) -conj(x(2)) -conj(x(3))];
        X_coded(7,:)  = [x(7) x(8) x(5) x(6) -conj(x(3)) -conj(x(4)) -conj(x(1)) -conj(x(2))];
        X_coded(8,:)  = [x(6) x(7) x(8) x(5) -conj(x(2)) -conj(x(3)) -conj(x(4)) -conj(x(1))];
        % Time slot 9–16 (conjugate phần tử)
        
        % Kênh truyền MIMO Rayleigh (8x8)
        Hr=(randn(frmLen,N,8)+ 1i*randn(frmLen,N,8))/sqrt(2);

        for n = 1:8
            %ma tran phat
            H = reshape(Hr(:,:,n),frmLen,8);
            %tin hieu tai anten thu
            r1(:,n)=awgn(sum(H.*X_coded(:,1),2)/sqrt(N),EbN0(idx));
            r2(:,n)=awgn(sum(H.*X_coded(:,2),2)/sqrt(N),EbN0(idx));
            r3(:,n)=awgn(sum(H.*X_coded(:,3),2)/sqrt(N),EbN0(idx));
            r4(:,n)=awgn(sum(H.*X_coded(:,4),2)/sqrt(N),EbN0(idx));
            r5(:,n)=awgn(sum(H.*X_coded(:,5),2)/sqrt(N),EbN0(idx));
            r6(:,n)=awgn(sum(H.*X_coded(:,6),2)/sqrt(N),EbN0(idx));
            r7(:,n)=awgn(sum(H.*X_coded(:,7),2)/sqrt(N),EbN0(idx));
            r8(:,n)=awgn(sum(H.*X_coded(:,8),2)/sqrt(N),EbN0(idx));
            
            % For z1 (estimating x1)
            z1_1(:,n) = r1(:,n).*conj(H(:,1)) + r2(:,n).*conj(H(:,2)) + r3(:,n).*conj(H(:,3)) + r4(:,n).*conj(H(:,4));
            z1_2(:,n) = conj(r5(:,n)).*H(:,5) + conj(r6(:,n)).*H(:,6) + conj(r7(:,n)).*H(:,7) + conj(r8(:,n)).*H(:,8);
            z1(:,n) = z1_1(:,n) - z1_2(:,n);
            
            % For z2 (estimating x2)
            z2_1(:,n) = r1(:,n).*conj(H(:,2)) + r2(:,n).*conj(H(:,1)) + r3(:,n).*conj(H(:,2)) + r4(:,n).*conj(H(:,3));
            z2_2(:,n) = conj(r5(:,n)).*H(:,8) + conj(r6(:,n)).*H(:,5) + conj(r7(:,n)).*H(:,6) + conj(r8(:,n)).*H(:,7) ;
            z2(:,n) = z2_1(:,n) - z2_2(:,n);
            
            % For z3 (estimating x3)
            z3_1(:,n) = r1(:,n).*conj(H(:,3)) + r2(:,n).*conj(H(:,4)) + r3(:,n).*conj(H(:,1)) + r4(:,n).*conj(H(:,2));
            z3_2(:,n) = conj(r5(:,n)).*H(:,7) + conj(r6(:,n)).*H(:,8) + conj(r7(:,n)).*H(:,5) + conj(r8(:,n)).*H(:,6);
            z3(:,n) = z3_1(:,n) - z3_2(:,n);
            
            % For z4 (estimating x4)
            z4_1(:,n) = r1(:,n).*conj(H(:,4)) + r2(:,n).*conj(H(:,3)) + r3(:,n).*conj(H(:,4)) + r4(:,n).*conj(H(:,1));
            z4_2(:,n) = conj(r5(:,n)).*H(:,6) + conj(r6(:,n)).*H(:,7) + conj(r7(:,n)).*H(:,8) + conj(r8(:,n)).*H(:,5);
            z4(:,n) = z4_1(:,n) - z4_2(:,n);
            
            % For z5 (estimating x5)
            z5_1(:,n) = r1(:,n).*conj(H(:,5)) + r2(:,n).*conj(H(:,6)) + r3(:,n).*conj(H(:,7)) + r4(:,n).*conj(H(:,8));
            z5_2(:,n) = conj(r5(:,n)).*H(:,1) + conj(r6(:,n)).*H(:,2) + conj(r7(:,n)).*H(:,3) + conj(r8(:,n)).*H(:,4);
            z5(:,n) = z5_1(:,n) - z5_2(:,n);
            
            % For z6 (estimating x6)
            z6_1(:,n) = r1(:,n).*conj(H(:,6)) + r2(:,n).*conj(H(:,5)) + r3(:,n).*conj(H(:,6)) + r4(:,n).*conj(H(:,7));
            z6_2(:,n) = conj(r5(:,n)).*H(:,4) + conj(r6(:,n)).*H(:,1) + conj(r7(:,n)).*H(:,2) + conj(r8(:,n)).*H(:,3);
            z6(:,n) = z6_1(:,n) - z6_2(:,n);
            
            % For z7 (estimating x7)
            z7_1(:,n) = r1(:,n).*conj(H(:,7)) + r2(:,n).*conj(H(:,8)) + r3(:,n).*conj(H(:,5)) + r4(:,n).*conj(H(:,6));
            z7_2(:,n) = conj(r5(:,n)).*H(:,3) + conj(r6(:,n)).*H(:,4) + conj(r7(:,n)).*H(:,1) + conj(r8(:,n)).*H(:,2);
            z7(:,n) = z7_1(:,n) - z7_2(:,n);
            
            % For z8 (estimating x8)
            z8_1(:,n) = r1(:,n).*conj(H(:,8)) + r2(:,n).*conj(H(:,7)) + r3(:,n).*conj(H(:,8)) + r4(:,n).*conj(H(:,5));
            z8_2(:,n) = conj(r5(:,n)).*H(:,2) + conj(r6(:,n)).*H(:,3) + conj(r7(:,n)).*H(:,4) + conj(r8(:,n)).*H(:,1);
            z8(:,n) = z8_1(:,n) - z8_2(:,n);

        end
        r(:,1) = sum(z1,2);
        r(:,2) = sum(z2,2);
        r(:,3) = sum(z3,2);
        r(:,4) = sum(z4,2);
        r(:,5) = sum(z5,2);
        r(:,6) = sum(z6,2);
        r(:,7) = sum(z7,2);
        r(:,8) = sum(z8,2);
        demod22 = pskdemod(r, P, 0);
        error22(packetIdx) = biterr(demod22,data);
    end
        BER22(idx) = sum(error22)/(numPackets*frmLen);
        fprintf('Eb/N0 = %d dB, BER = = %e\n', EbN0(idx), BER22(idx));
        semilogy(EbN0(1:idx),BER22(1:idx),'r*-')
end
% Vẽ đồ thị
h = figure;
semilogy(EbN0, BER22, 'r*', 'LineWidth', 1.5, 'DisplayName', 'Tx=8,Rx=1 MIMO 8x8');
grid on;
set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-5 1]);
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('Hệ thống MIM0 8x8');
legend('show');
set(h, 'color', [1 1 1]);
set(h, 'NumberTitle', 'off');
set(h, 'renderer', 'zbuffer');