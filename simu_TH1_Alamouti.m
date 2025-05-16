function simu_TH1_Alamouti(frmLen,numPackets,EbN0dB)
    numTx = 2;             
    numRx1 = 1;            
    numRx2 = 2;           
    M = 4;                 
    for idx = 1:length(EbN0dB)
        for pktIdx = 1:numPackets
            data = randi([0 M-1], frmLen, 2);  
            modData = pskmod(data, M, 0);      
            dataSISO = data(:, 1);
         
            %SISO
            xSISO = pskmod(dataSISO, M, 0);            
            hSISO = (randn(frmLen, 1) + 1i*randn(frmLen, 1)) / sqrt(2);
            

            SNR_dB = EbN0dB(idx) + 10*log10(log2(M));
            rSISO = hSISO .* xSISO;
            rSISO = awgn(rSISO, SNR_dB,'measured'); 
            
            xHatSISO = conj(hSISO) .* rSISO ./ (abs(hSISO).^2);
            decodedDataSISO = pskdemod(xHatSISO, M, 0);
            
            bitErrSISO(pktIdx) = biterr(decodedDataSISO, dataSISO);
    
            %Alamouti 2x1
            x1 = modData(:,1);
            x2 = modData(:,2);
            
            H = (randn(frmLen, numTx, numRx1) + 1i*randn(frmLen, numTx, numRx1)) / sqrt(2);
            h1 = H(:,1,1);
            h2 = H(:,2,1);
    
            s1 = h1 .* x1 + h2 .* x2;
            s2 = -h1 .* conj(x2) + h2 .* conj(x1);
    
            r1 = awgn(s1, SNR_dB,'measured');
            r2 = awgn(s2, SNR_dB,'measured');
    
            x1Hat = conj(h1).*r1 + h2.*conj(r2);
            x2Hat = conj(h2).*r1 - h1.*conj(r2);
    
            norm_factor = sum(sum(abs(H).^2, 3), 2);
            decodedSymbolsAla21 = [x1Hat x2Hat]./norm_factor;
            decodedDataAla21 = pskdemod(decodedSymbolsAla21, M, 0);
            
            bitErrAla21(pktIdx) = biterr(decodedDataAla21, data);
    

            %Alamouti 2x2
            s1 = modData(:,1);
            s2 = modData(:,2);
    
            tx1 = [s1 s2];
            tx2 = [-conj(s2) conj(s1)];
    
            H = (randn(frmLen, numTx, numRx2) + 1i*randn(frmLen, numTx, numRx2)) / sqrt(2);
            h11 = H(:,1,1); h12 = H(:,2,1);
            h21 = H(:,1,2); h22 = H(:,2,2);
    
            r11 = awgn(h11 .* s1 + h12 .* s2,  SNR_dB,'measured');
            r12 = awgn(h11 .* (-conj(s2)) + h12 .* conj(s1), SNR_dB,'measured');
            r21 = awgn(h21 .* s1 + h22 .* s2, SNR_dB,'measured');
            r22 = awgn(h21 .* (-conj(s2)) + h22 .* conj(s1), SNR_dB,'measured');
   
            z1 = conj(h11).*r11 + h12.*conj(r12) + conj(h21).*r21 + h22.*conj(r22);
            z2 = conj(h12).*r11 - h11.*conj(r12) + conj(h22).*r21 - h21.*conj(r22);
   
            z(:,1) = sum(z1,2);
            z(:,2) = sum(z2,2);
    
            norm_factor = sum(sum(abs(H).^2, 3), 2);
            z=z./norm_factor;
            decodedDataAla22 = pskdemod(z, M, 0);

            bitErrAla22(pktIdx) = biterr(decodedDataAla22, data);
        end
    

        berSISO(idx) = sum(bitErrSISO) / (numPackets * frmLen * log2(M));
        berAlamouti21(idx) = sum(bitErrAla21) / (numPackets * frmLen * 2 * log2(M));
        berAlamouti22(idx) = sum(bitErrAla22) / (numPackets * frmLen * 2 * log2(M));
        
        fprintf('Eb/N0 = %d dB, BER Alamouti 2x1 = %.3e, BER Alamouti 2x2 = %.3e, BER SISO = %.3e\n', ...
                EbN0dB(idx), berAlamouti21(idx), berAlamouti22(idx), berSISO(idx));
    end
    
    figure;
    semilogy(EbN0dB, berSISO, 'g^-', 'LineWidth', 1.5, 'DisplayName', 'SISO 1x1'); hold on;
    semilogy(EbN0dB, berAlamouti21, 'r*-', 'LineWidth', 1.5, 'DisplayName', 'Alamouti 2x1');
    semilogy(EbN0dB, berAlamouti22, 'bs-', 'LineWidth', 1.5, 'DisplayName', 'Alamouti 2x2');
    grid on;
    
    set(gca, 'yscale', 'log', 'xlim', [EbN0dB(1) EbN0dB(end)], 'ylim', [1e-8 1]);
    xlabel('Eb/N0 (dB)');
    ylabel('Bit Error Rate (BER)');
    title('BER vs SNR for Alamouti with QPSK');
    legend('show');
    set(gcf, 'color', [1 1 1]);
end
