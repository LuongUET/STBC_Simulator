%Mo phong TH1
frmLen     = 1e5;
numPackets = 100;
EbN0       = 0:2:20;
simu_TH1_Alamouti(frmLen,numPackets,EbN0);

%Mo phong TH2
frmLen     = 1e5;
numPackets = 100;
EbN0       = 0:2:20;
%simu_TH2_STBC3x4(frmLen,numPackets,EbN0);

%Mo phong TH3
frmLen     = 1e5;
numPackets = 100;
EbN0       = 0:2:20;
%simu_TH3_STBC4x4(frmLen,numPackets,EbN0);

%Mo phong TH4
frmLen     = 1e5; 
numPackets = 100;     
EbN0       = 0:2:20; 
%simu_TH4_STBC4x4(frmLen,numPackets,EbN0);