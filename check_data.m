clear;
close all;
raw_data = load('raw_data.txt');
raw_data = raw_data-1800;

filter_data = load('filter_data.txt');
filter_data = filter_data-1800;
sd_rd=std(raw_data);
sd_rdf=std(filter_data);

[l,N]=size(raw_data);  %N为数据长度
%raw_data= raw_data(99:N);
figure(1);
plot(raw_data);
hold on;
plot(filter_data);
y_data = load('m_ydata.txt');
x_data = load('m_xdata.txt');
hold off;
 
figure;
plot(raw_data);
hold on;
%plot(x_data,y_data,'-+');hold off;
plot(x_data,y_data,'o');hold off;
% 显示滤波后的极值
yf_data = load('m_ydata_f.txt');
xf_data = load('m_xdata_f.txt');
figure;
plot(filter_data);
hold on;
%plot(x_data,y_data,'-+');hold off;
plot(xf_data,yf_data,'o');hold off;
%peaks = AMPD_P(filter_data,2,1);

 fs=250000;         %采样率 % Sampling Frequency 
%NFFT=2^nextpow2(N);  
NFFT=N;
n=0:N-1;
t=n/fs;
f=(0:N-1)*fs/N;

y_raw=fft(raw_data,NFFT);
mag1=abs(y_raw);



y_filter=fft(filter_data,NFFT);
mag2=abs(y_filter);

figure(4);
plot(f(1:N/2),mag1(1:N/2)*2/N ,f(1:N/2),mag2(1:N/2)*2/N );  
 
title('FFT 频谱');
xlabel('频率/HZ');ylabel('振幅');
axis([000 60000 0 100]); 


% figure(2);
% subplot(211);
% plot(f(1:N/2),mag1(1:N/2)*2/N);
% title('FFT 频谱');
% xlabel('频率/HZ');ylabel('振幅');
% axis([0 1000 0 3000]); 
% subplot(212);
% plot(f(1:N/2),mag2(1:N/2)*2/N);
% title('FFT 频谱');
% xlabel('频率/HZ');ylabel('振幅');
% axis([0 1000 0 3000]); 
