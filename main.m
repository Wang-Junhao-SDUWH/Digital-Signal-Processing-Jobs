clc;clear;
%%读入源文件
[y,FS] = audioread('D:\University\作业\数字信号处理\音频\spring.m4a');
%合并声道
if size(y,2)==2
    y = (y(:,1)+y(:,2))/2;
end
y = y/max(abs(y)); %归一化
T = 1/FS;
t = (0:length(y)-1)*T;

%%语音信号预处理
%1.消除直流分量
y = y - mean(y);
%2.滤波▲
%切比雪夫Ⅱ型高通滤波器
wp=2*50/FS;
ws=2*55/FS;
rp =3;rs = 20;
[N,Wn] = cheb2ord(wp,ws,rp,rs);
[b,a] = cheby2(N,rs,Wn);
y = filter(b,a,y);
y = y/max(abs(y));

subplot 211;
plot(t,y,'b');
title("预处理后的声音");
xlabel("时间");ylabel("音量");
sound(y,FS);

output = pitch_alter(y,FS);

subplot 212;
plot(t,output,'b');
title("变调后的声音");
ylabel("音量");xlabel("时间");
sound(output,FS);
wlen = FS*25/1000;       %以25ms为一帧
inc = 2*wlen/5;          %帧移