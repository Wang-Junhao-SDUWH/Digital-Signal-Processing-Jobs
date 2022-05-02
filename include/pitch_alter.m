function output = pitch_alter(y,FS)
    %%提取信息
    wlen = FS*25/1000;       %以25ms为一帧
    inc = 2*wlen/5;          %帧移
    overlap = wlen - inc;   %重叠长度
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
    p = 12; %设预测阶次
    frameTime = frame2time(fn,wlen,inc,FS); %每帧对应的时间刻度
    rate = input('请输入变声力度\n (0,1)->男化 \n （1，+∞)->女化\n');
    %rate = str2num(in);
    for i=1:fn          %求取每帧的预测系数和增益
        u=Y(:,i);
        [ar,g]=lpc(u,p);
        Ar_coe(:,i) = ar;
        Gain(i) = g;
    end
    %基音检测
    [cbase1,fbase1,~,SF,~,~,~,~,~]=Ext_F0ztms(y,FS,wlen,inc,T1,r2,minL,minlo,ThrC,0);
    if rate>1
        sign = -1;
    else
        sign = 1;
    end
    lmin = floor(FS/450);           %基音周期的最小值
    lmax = floor(FS/60);            %基音周期的最大值
    deltaOMG = sign*100*2*pi/FS;    %根植顺时针或逆时针旋转量
    cbase2 = cbase1/rate;          %增减后的基音周期
    fbase2 = fbase1*rate;            %增减后的基音频率

    m = 0;    %初始化
    zint = zeros(p,1);
    for i=1:fn
        a = Ar_coe(:,i);      %取得本帧的AR系数
        sigma = Gain(i);  %取得本帧的增益系数
        sigma = sqrt(sigma);
        if SF(i)==0           %无话帧
            excitation = randn(wlen,1); %产生白噪声
            [synt_frame,zint]=filter(sigma,a,excitation,zint);
        else      %有话帧
            tempT = floor(cbase2(i)); %把周期值变为整数
            if tempT<lmin,tempT = lmin;end
            if tempT>lmax,tempT = lmax;end
            ft=roots(a);
            ft1=ft;
            %增加共振峰频率，实轴上方的根顺时针转，下方的根逆时针转，求出新的根植
            for k=1:p
                if imag(ft(k)>0)
                    ft1(k) = ft(k)*exp(1i*deltaOMG);
                elseif imag(ft(k))<0
                    ft1(k) = ft(k)*exp(-1i*deltaOMG);
                end
            end
            ai = poly(ft1); %由新的根植重新组成预测系数
            exc_syn1 = zeros(wlen+m,1);   %初始化脉冲发生区
            exc_syn1(mod(1:m+wlen,tempT)==0)=1;  %在基音周期的位置产生脉冲，幅值为1
            exc_syn2 = exc_syn1(m+1:m+inc); %计算帧移inc区间内的脉冲个数
            index = find(exc_syn2==1);
            excitation = exc_syn1(m+1:m+wlen);  %这一帧的激励脉冲源
            if isempty(index)
                m = m+inc;
            else
                eal = length(index);    %计算脉冲个数
                m = inc-index(eal);   %计算下一帧的前导零点
            end
            gain = sigma/sqrt(1/tempT);    %增益
            [synt_frame,zint] = filter(gain,ai,excitation,zint);%用激励脉冲合成语音
        end
        if i==1
            output = synt_frame;
        else
            M = length(output); %重叠部分的处理
            output = [output(1:M-overlap);output(M-overlap+1:M).*...
                temp1+synt_frame(1:overlap).*temp2;synt_frame(overlap+1:wlen)];
        end
    end
    output(isnan(output))=0;
    lo = length(output);    %把输出output延长至与输入信号xx等长
    if lo<ly
        tempo = [output;zeros(ly-lo,1)];
    else
        tempo = output(1:ly);
    end
    output = real(tempo/max(abs(tempo)));
end

