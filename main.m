clc;clear;
%%����Դ�ļ�
[y,FS] = audioread('D:\University\��ҵ\�����źŴ���\��Ƶ\spring.m4a');
%�ϲ�����
if size(y,2)==2
    y = (y(:,1)+y(:,2))/2;
end
y = y/max(abs(y)); %��һ��
T = 1/FS;
t = (0:length(y)-1)*T;

%%�����ź�Ԥ����
%1.����ֱ������
y = y - mean(y);
%2.�˲���
%�б�ѩ����͸�ͨ�˲���
wp=2*50/FS;
ws=2*55/FS;
rp =3;rs = 20;
[N,Wn] = cheb2ord(wp,ws,rp,rs);
[b,a] = cheby2(N,rs,Wn);
y = filter(b,a,y);
y = y/max(abs(y));

subplot 211;
plot(t,y,'b');
title("Ԥ����������");
xlabel("ʱ��");ylabel("����");
sound(y,FS);

output = pitch_alter(y,FS);

subplot 212;
plot(t,output,'b');
title("����������");
ylabel("����");xlabel("ʱ��");
sound(output,FS);
wlen = FS*25/1000;       %��25msΪһ֡
inc = 2*wlen/5;          %֡��