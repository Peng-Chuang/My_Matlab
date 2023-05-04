% Lf_road_num = 1;  %低频通道
% 
% Lf_data_1 = dataWave(110324:110324+20063);
% Lfdata=Lf_data_1(:,Lf_road_num);  %第几通道
% fir_data_out=fir_calc(0,Lfdata);
% [window_data_out]=window_fixed_point(1,fir_data_out);
% [fft_data_out]=fft_fixed_point(window_data_out);
% fft_data_out_Lf=fft_data_out(2:401);
% 
% figure(1)
% plot(fft_data_out_Lf);
% hold on
% %plot(Lf_feature((Lf_road_num-1)*400+1:400*Lf_road_num));
% % plot(Lf_feature(401:800));  %低频计算值  量纲  3e4   一毫安
% title("低频");
% grid on;

Hf_data_1=dataWave(21436+270+69+10000:21436+270+69+30064);%20064个点，20ms采样+2msFFT计算
Hf_data_1=Hf_data_1(500:end);
dianshu = 20064;
time = (1/50)/2;
step= time*10^6;
Hz=700;
T=1/(2*Hz);
detaT = 0.000001;
DCloudian=[];ACloudian=[];
numb =dianshu-step+1;

    W=4000;
    Xsequence1=Hf_data_1(W:W+1500-1);
    XMM2=Xsequence1(1:1:end);
%     XMM2 = reshape(XMM,100,[]);
%     XMM2 = mean(XX,1);
%     XMM2=XMM(1:10:end);
 
    
    nn=0:length(XMM2)-1;
    H = [sin(2*50*pi*detaT*nn)',cos(2*50*pi*detaT*nn)',ones(length(nn),1)];
    X =inv(H'*H)*H'* XMM2;
    A=sqrt(X(1)^2+X(2)^2);
    B=X(3);
    DCloudian=[DCloudian;B];
    ACloudian =[ACloudian;A];



% nn=0:length(Hf_data_1)-1;
% 
% func_sin = @(a,t) a(1)*sin(2*pi*50*T*t/50+a(2)) + a(3);
% 
% Aa= lsqcurvefit( func_sin, [5000 800 2000], nn, Hf_data_1');

