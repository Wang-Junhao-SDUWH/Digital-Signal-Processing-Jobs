function output = pitch_alter(y,FS)
    %%��ȡ��Ϣ
    wlen = FS*25/1000;       %��25msΪһ֡
    inc = 2*wlen/5;          %֡��
    overlap = wlen - inc;   %�ص�����
    temp1 = (0:overlap-1)'/overlap;
    temp2 = (overlap-1:-1:0)'/ overlap;
    n2 = 1:wlen/2+1;
    Y = enframe(y,wlen,inc)';
    ly = length(y);
    fn = size(Y,2);
    T1 = 0.1;r2 = 0.5;
    minL = 10;
    minlo = 5;
    ThrC = [10 15];
    p = 12; %��Ԥ��״�
    frameTime = frame2time(fn,wlen,inc,FS); %ÿ֡��Ӧ��ʱ��̶�
    rate = input('�������������\n (0,1)->�л� \n ��1��+��)->Ů��\n');
    %rate = str2num(in);
    for i=1:fn          %��ȡÿ֡��Ԥ��ϵ��������
        u=Y(:,i);
        [ar,g]=lpc(u,p);
        Ar_coe(:,i) = ar;
        Gain(i) = g;
    end
    %�������
    [cbase1,fbase1,~,SF,~,~,~,~,~]=Ext_F0ztms(y,FS,wlen,inc,T1,r2,minL,minlo,ThrC,0);
    if rate>1
        sign = -1;
    else
        sign = 1;
    end
    lmin = floor(FS/450);           %�������ڵ���Сֵ
    lmax = floor(FS/60);            %�������ڵ����ֵ
    deltaOMG = sign*100*2*pi/FS;    %��ֲ˳ʱ�����ʱ����ת��
    cbase2 = cbase1/rate;          %������Ļ�������
    fbase2 = fbase1*rate;            %������Ļ���Ƶ��

    m = 0;    %��ʼ��
    zint = zeros(p,1);
    for i=1:fn
        a = Ar_coe(:,i);      %ȡ�ñ�֡��ARϵ��
        sigma = Gain(i);  %ȡ�ñ�֡������ϵ��
        sigma = sqrt(sigma);
        if SF(i)==0           %�޻�֡
            excitation = randn(wlen,1); %����������
            [synt_frame,zint]=filter(sigma,a,excitation,zint);
        else      %�л�֡
            tempT = floor(cbase2(i)); %������ֵ��Ϊ����
            if tempT<lmin,tempT = lmin;end
            if tempT>lmax,tempT = lmax;end
            ft=roots(a);
            ft1=ft;
            %���ӹ����Ƶ�ʣ�ʵ���Ϸ��ĸ�˳ʱ��ת���·��ĸ���ʱ��ת������µĸ�ֲ
            for k=1:p
                if imag(ft(k)>0)
                    ft1(k) = ft(k)*exp(1i*deltaOMG);
                elseif imag(ft(k))<0
                    ft1(k) = ft(k)*exp(-1i*deltaOMG);
                end
            end
            ai = poly(ft1); %���µĸ�ֲ�������Ԥ��ϵ��
            exc_syn1 = zeros(wlen+m,1);   %��ʼ�����巢����
            exc_syn1(mod(1:m+wlen,tempT)==0)=1;  %�ڻ������ڵ�λ�ò������壬��ֵΪ1
            exc_syn2 = exc_syn1(m+1:m+inc); %����֡��inc�����ڵ��������
            index = find(exc_syn2==1);
            excitation = exc_syn1(m+1:m+wlen);  %��һ֡�ļ�������Դ
            if isempty(index)
                m = m+inc;
            else
                eal = length(index);    %�����������
                m = inc-index(eal);   %������һ֡��ǰ�����
            end
            gain = sigma/sqrt(1/tempT);    %����
            [synt_frame,zint] = filter(gain,ai,excitation,zint);%�ü�������ϳ�����
        end
        if i==1
            output = synt_frame;
        else
            M = length(output); %�ص����ֵĴ���
            output = [output(1:M-overlap);output(M-overlap+1:M).*...
                temp1+synt_frame(1:overlap).*temp2;synt_frame(overlap+1:wlen)];
        end
    end
    output(isnan(output))=0;
    lo = length(output);    %�����output�ӳ����������ź�xx�ȳ�
    if lo<ly
        tempo = [output;zeros(ly-lo,1)];
    else
        tempo = output(1:ly);
    end
    output = real(tempo/max(abs(tempo)));
end

