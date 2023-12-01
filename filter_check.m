clear;
close all;
raw_data = load('raw_data.txt');
raw_data = raw_data-1800;

filter_data = load('filter_data.txt');
filter_data = filter_data-1800;
 
sd_rd=std(raw_data,'omitnan');
sd_rdf1=std(filter_data);

[l,N]=size(raw_data);  %N为数据长度

figure;
plot(raw_data);
hold on;
plot(filter_data);
hold off;

order = 3;
framelen =5;
filter_d2 = sgolayfilt(raw_data,order,framelen);
sd_rdf2=std(filter_d2,'omitnan');

figure;
plot(raw_data);
hold on;
plot(filter_d2);
hold off;
figure;
plot(filter_data);
hold on;
plot(filter_d2);
hold off;


 fs=250000;         %采样率 % Sampling Frequency 
NFFT=2^nextpow2(N);  
%NFFT=N;
n=0:N-1;
t=n/fs;
f=(0:N-1)*fs/N;

y_raw=fft(raw_data,NFFT);
mag1=abs(y_raw);



y_filter=fft(filter_data,NFFT);
mag2=abs(y_filter);

figure;
plot(f(1:N/2),mag1(1:N/2)*2/N ,f(1:N/2),mag2(1:N/2)*2/N );  
 
title('FFT 频谱');
xlabel('频率/HZ');ylabel('振幅');
axis([000 60000 0 100]); 
