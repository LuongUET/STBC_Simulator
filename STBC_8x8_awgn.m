clc; clear;

N = 8; % Số anten phát
M = 8; % Số anten thu
P = 4; % QPSK
EbN0_dB = 0:2:20;
frmLen = 8; % Số frame để tính BER
BER = zeros(1, length(EbN0_dB));

for idx = 1:length(EbN0_dB)
    errors = 0; total = 0;
    snr = EbN0_dB(idx);
    
    for pkt = 1:frmLen
        %Tạo dữ liệu và điều chế QPSK
        data = randi([0 P-1], frmLen, 8); % 8 symbol
        x = pskmod(data, P, pi/4);  % QPSK
        
        %Mã hóa STBC theo G8
        X_coded = zeros(16, 8); % 16 time slot × 8 anten
        % Time slot 1–8
        X_coded(1,:)  = [x(1)  x(2)  x(3)  x(4)  x(5)  x(6)  x(7)  x(8)];
        X_coded(2,:)  = [-x(2) x(1)  x(4) -x(3)  x(6) -x(5) -x(8)  x(7)];
        X_coded(3,:)  = [-x(3) -x(4) x(1)  x(2)  x(7)  x(8) -x(5) -x(6)];
        X_coded(4,:)  = [-x(4) x(3) -x(2) x(1)  x(8) -x(7)  x(6) -x(5)];
        X_coded(5,:)  = [-x(5) -x(6) -x(7) -x(8) x(1)  x(2)  x(3)  x(4)];
        X_coded(6,:)  = [-x(6) x(5) -x(8) x(7) -x(2) x(1) -x(4)  x(3)];
        X_coded(7,:)  = [-x(7) x(8) x(5) -x(6) -x(3)  x(4) x(1) -x(2)];
        X_coded(8,:)  = [-x(8) -x(7) x(6) x(5) -x(4) -x(3)  x(2)  x(1)];
        % Time slot 9–16 (conjugate phần tử)
        for i = 1:8
            X_coded(8+i,:) = conj(X_coded(i,:));
        end
        
        % Kênh truyền MIMO Rayleigh (8x8)
        Hr=(randn(frmLen*2,N,8)+ 1i*randn(frmLen*2,N,8))/sqrt(2);

        for n = 1:8
            %ma tran phat
            H = reshape(Hr(:,:,n),frmLen*2,N);
            %tin hieu tai anten thu
            r1(:,n)=awgn(sum(H.*X_coded(:,1),2)/sqrt(N),EbN0(idx));
            r2(:,n)=awgn(sum(H.*X_coded(:,2),2)/sqrt(N),EbN0(idx));
            r3(:,n)=awgn(sum(H.*X_coded(:,3),2)/sqrt(N),EbN0(idx));
            r4(:,n)=awgn(sum(H.*X_coded(:,4),2)/sqrt(N),EbN0(idx));
            r5(:,n)=awgn(sum(H.*X_coded(:,5),2)/sqrt(N),EbN0(idx));
            r6(:,n)=awgn(sum(H.*X_coded(:,6),2)/sqrt(N),EbN0(idx));
            r7(:,n)=awgn(sum(H.*X_coded(:,7),2)/sqrt(N),EbN0(idx));
            r8(:,n)=awgn(sum(H.*X_coded(:,8),2)/sqrt(N),EbN0(idx));
            r9(:,n)=awgn(sum(H.*X_coded(:,9),2)/sqrt(N),EbN0(idx));
            r10(:,n)=awgn(sum(H.*X_coded(:,10),2)/sqrt(N),EbN0(idx));
            r11(:,n)=awgn(sum(H.*X_coded(:,11),2)/sqrt(N),EbN0(idx));
            r12(:,n)=awgn(sum(H.*X_coded(:,12),2)/sqrt(N),EbN0(idx));
            r13(:,n)=awgn(sum(H.*X_coded(:,13),2)/sqrt(N),EbN0(idx));
            r14(:,n)=awgn(sum(H.*X_coded(:,14),2)/sqrt(N),EbN0(idx));
            r15(:,n)=awgn(sum(H.*X_coded(:,15),2)/sqrt(N),EbN0(idx));
            r16(:,n)=awgn(sum(H.*X_coded(:,16),2)/sqrt(N),EbN0(idx));
            
            % For z1 (estimating x1)
            z1_1(:,n) = r1(:,n).*conj(H(:,1)) + r2(:,n).*conj(H(:,2)) + r3(:,n).*conj(H(:,3)) + r4(:,n).*conj(H(:,4));
            z1_2(:,n) = conj(r5(:,n)).*H(:,5) + conj(r6(:,n)).*H(:,6) + conj(r7(:,n)).*H(:,7) + conj(r8(:,n)).*H(:,8);
            z1(:,n) = z1_1(:,n) - z1_2(:,n);
            
            % For z2 (estimating x2)
            z2_1(:,n) = r2(:,n).*conj(H(:,1)) + r3(:,n).*conj(H(:,2)) + r4(:,n).*conj(H(:,3)) + r1(:,n).*conj(H(:,4));
            z2_2(:,n) = conj(r6(:,n)).*H(:,5) + conj(r7(:,n)).*H(:,6) + conj(r8(:,n)).*H(:,7) + conj(r5(:,n)).*H(:,8);
            z2(:,n) = z2_1(:,n) - z2_2(:,n);
            
            % For z3 (estimating x3)
            z3_1(:,n) = r3(:,n).*conj(H(:,1)) + r4(:,n).*conj(H(:,2)) + r1(:,n).*conj(H(:,3)) + r2(:,n).*conj(H(:,4));
            z3_2(:,n) = conj(r7(:,n)).*H(:,5) + conj(r8(:,n)).*H(:,6) + conj(r5(:,n)).*H(:,7) + conj(r6(:,n)).*H(:,8);
            z3(:,n) = z3_1(:,n) - z3_2(:,n);
            
            % For z4 (estimating x4)
            z4_1(:,n) = r4(:,n).*conj(H(:,1)) + r1(:,n).*conj(H(:,2)) + r2(:,n).*conj(H(:,3)) + r3(:,n).*conj(H(:,4));
            z4_2(:,n) = conj(r8(:,n)).*H(:,5) + conj(r5(:,n)).*H(:,6) + conj(r6(:,n)).*H(:,7) + conj(r7(:,n)).*H(:,8);
            z4(:,n) = z4_1(:,n) - z4_2(:,n);
            
            % For z5 (estimating x5)
            z5_1(:,n) = r5(:,n).*conj(H(:,1)) + r6(:,n).*conj(H(:,2)) + r7(:,n).*conj(H(:,3)) + r8(:,n).*conj(H(:,4));
            z5_2(:,n) = conj(r1(:,n)).*H(:,5) + conj(r2(:,n)).*H(:,6) + conj(r3(:,n)).*H(:,7) + conj(r4(:,n)).*H(:,8);
            z5(:,n) = -z5_1(:,n) + z5_2(:,n);
            
            % For z6 (estimating x6)
            z6_1(:,n) = r6(:,n).*conj(H(:,1)) + r7(:,n).*conj(H(:,2)) + r8(:,n).*conj(H(:,3)) + r5(:,n).*conj(H(:,4));
            z6_2(:,n) = conj(r2(:,n)).*H(:,5) + conj(r3(:,n)).*H(:,6) + conj(r4(:,n)).*H(:,7) + conj(r1(:,n)).*H(:,8);
            z6(:,n) = -z6_1(:,n) + z6_2(:,n);
            
            % For z7 (estimating x7)
            z7_1(:,n) = r7(:,n).*conj(H(:,1)) + r8(:,n).*conj(H(:,2)) + r5(:,n).*conj(H(:,3)) + r6(:,n).*conj(H(:,4));
            z7_2(:,n) = conj(r3(:,n)).*H(:,5) + conj(r4(:,n)).*H(:,6) + conj(r1(:,n)).*H(:,7) + conj(r2(:,n)).*H(:,8);
            z7(:,n) = -z7_1(:,n) + z7_2(:,n);
            
            % For z8 (estimating x8)
            z8_1(:,n) = r8(:,n).*conj(H(:,1)) + r5(:,n).*conj(H(:,2)) + r6(:,n).*conj(H(:,3)) + r7(:,n).*conj(H(:,4));
            z8_2(:,n) = conj(r4(:,n)).*H(:,5) + conj(r1(:,n)).*H(:,6) + conj(r2(:,n)).*H(:,7) + conj(r3(:,n)).*H(:,8);
            z8(:,n) = -z8_1(:,n) + z8_2(:,n);

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
semilogy(EbN0, BER22, 'r*-', 'LineWidth', 1.5, 'DisplayName', 'MIMO 4x4');
grid on;
set(gca, 'yscale', 'log', 'xlim', [EbN0(1), EbN0(end)], 'ylim', [1e-5 1]);
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('Hệ thống MIM0 4x4');
legend('show');
set(h, 'color', [1 1 1]);
set(h, 'NumberTitle', 'off');
set(h, 'renderer', 'zbuffer');