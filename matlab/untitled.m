clear all
close all
set(0,'defaultfigurecolor','w')

l1=[9.8,10.05,5.135,5.135,0.09,0.09,4.97,4.97,9.25,9.25,10.455,10.455,10.455,6.165,6.165,6.165;9.825,8.59,3.66,3.66,1.58,1.58,6.475,6.475,10.765,10.765,7.50,7.50,7.50,3.20,3.20,3.20]
l2=l1+0.02;
fstart=[2,4,9,9,14,14,19,19,24,24,31,31,31,37,37,37]*1e6;
fend=[8,10,18,18,26,26,34,34,42,42,48,48,48,55,55,55]*1e6

%%
% 

fp=152525;

stp=500825;

str=0.05*stp+fp; sto=0.45*stp+fp; 

str1=0.55*stp+fp; sto1=0.95*stp+fp;



tag=[1 2 3 5 7 9 11 14]


for m=1:length(tag)

j=tag(m);  


% L1=3.98; L2=4.00%调焦系统，铁板目标
mar_zero=[];
martmp=[];nor_smo=[];nortmp=[];tmp=[];tmp_mul=[];mar_peak=[];
for i=1:2
if(j<=10)
if(i==1)
xnfs0=load(['H:\tsinghua\laser lidar\202601\20260902-多短臂长差辅助干涉光路降频\5-实际目标测试\1023_00' num2str(j-1) '_ch2.mat']);
else
xnfs0=load(['H:\tsinghua\laser lidar\202601\20260902-多短臂长差辅助干涉光路降频\5-实际目标测试\1023_00' num2str(j-1) '_ch3.mat']);
end
xns0=load(['H:\tsinghua\laser lidar\202601\20260902-多短臂长差辅助干涉光路降频\5-实际目标测试\1023_00' num2str(j-1) '_ch4.mat']);
xnfs0=xnfs0.data;
xns0=xns0.data;
else
if(i==1)
xnfs0=load(['H:\tsinghua\laser lidar\202601\20260902-多短臂长差辅助干涉光路降频\5-实际目标测试\1023_0' num2str(j-1) '_ch2.mat']);
else
xnfs0=load(['H:\tsinghua\laser lidar\202601\20260902-多短臂长差辅助干涉光路降频\5-实际目标测试\1023_0' num2str(j-1) '_ch3.mat']);
end

xns0=load(['H:\tsinghua\laser lidar\202601\20260902-多短臂长差辅助干涉光路降频\5-实际目标测试\1023_0' num2str(j-1) '_ch4.mat']);
xnfs0=xnfs0.data;
xns0=xns0.data;
end

xnfs0=xnfs0(:,1);
figure, plot(xnfs0());
xns0=xns0(:,1);
L1=l1(i,j);
L2=l2(i,j);
 if(i==1)
        Laux1= 24.964926;
        
       


 else
        Laux2= 23.539728;
        


 end








xnfs0=(xnfs0-mean(xnfs0));
xns0=(xns0-mean(xns0));

%%
if(i==1)
    shift=0;
else
    shift=0;
end
for k=1:10
xnfs=xnfs0(str+stp*(k-1):sto+stp*(k-1));
xnfs1=xnfs0(str1+stp*(k-1):sto1+stp*(k-1));


xns=xns0(str+stp*(k-1)-shift:sto+stp*(k-1)-shift);
xns1=xns0(str1+stp*(k-1)-shift:sto1+stp*(k-1)-shift);



fs=250e6;
nyq = fs / 2;
if(i==1)
bandpassHz = [2e6, 14e6];
else
bandpassHz = [2e6, 14e6];
end
[b, a] = butter(4, bandpassHz / nyq, 'bandpass');
xnfs = filtfilt(b, a, xnfs);
xnfs=(xnfs-mean(xnfs));
xnfs=xnfs./abs(hilbert(xnfs));
xnfs1 = filtfilt(b, a, xnfs1);
xnfs1=(xnfs1-mean(xnfs1));
xnfs1=xnfs1./abs(hilbert(xnfs1));


xnfs=xnfs(0.01*length(xnfs):0.99*length(xnfs));
xnfs1=xnfs1(0.01*length(xnfs1):0.99*length(xnfs1));



% 
bandpassHz2 = [fstart(1,j), fend(1,j)];
[b2, a2] = butter(4, bandpassHz2 / nyq, 'bandpass');
xns = filtfilt(b2, a2, xns);
xns=(xns-mean(xns));
% xns=xns./abs(hilbert(xns));
xns1 = filtfilt(b2, a2, xns1);
xns1=(xns1-mean(xns1));
% xns1=xns1./abs(hilbert(xns1));
xns=xns(0.01*length(xns):0.99*length(xns));
xns1=xns1(0.01*length(xns1):0.99*length(xns1));

% env_safe = abs(hilbert(xns));
% safe_threshold = 0.01 * max(env_safe); 
% env_safe(env_safe < safe_threshold) = safe_threshold;
% xns = xns ./ env_safe;
% 
% env_safe1 = abs(hilbert(xns1));
% safe_threshold1 = 0.01 * max(env_safe1); 
% env_safe1(env_safe1 < safe_threshold1) = safe_threshold1;
% xns1 = xns1 ./ env_safe1;


[xns_cal_peak,mar_peak]=ReSample20260908(xnfs,xns,1);%！
[xns_cal_peak1,mar_peak1]=ReSample20260908(xnfs1,xns1,1);%！


% figure, hold on
% plot(xnfs)
% ylim([-2 2])
% ttmp=1:length(xnfs);
% plot(ttmp(mar_peak),xnfs(mar_peak),'ro')
% figure, hold on
% plot(xnfs1)
% ylim([-2 2])
% ttmp=1:length(xnfs1);
% plot(ttmp(mar_peak1),xnfs1(mar_peak1),'ro')



if(i==1)
[freqout,specout0,lcmax,pcmax]=CZT(xns_cal_peak,L1,L2,200000,Laux1);
pktmp=peakget(freqout,L1,L2,20*log10(abs(specout0)));
RS1(k,j)=pktmp; %原始信号 ,用作对比

[freqout1,specout1,lcmax,pcmax]=CZT(xns_cal_peak1,L1,L2,200000,Laux1);
pktmp=peakget(freqout1,L1,L2,20*log10(abs(specout1)));
RS2(k,j)=pktmp; %原始信号 ,用作对比
RS(k,j)=RS1(k,j)/2+RS2(k,j)/2;
span(j)=(length(xns_cal_peak)+length(xns_cal_peak1))/2*3e8/2/Laux1;
else
[freqout,specout0,lcmax,pcmax]=CZT(xns_cal_peak,L1,L2,200000,Laux2);
pktmp=peakget(freqout,L1,L2,20*log10(abs(specout0)));
RS1(k,j)=pktmp; %原始信号 ,用作对比

[freqout1,specout1,lcmax,pcmax]=CZT(xns_cal_peak1,L1,L2,200000,Laux2);
pktmp=peakget(freqout1,L1,L2,20*log10(abs(specout1)));
RS2(k,j)=pktmp; %原始信号 ,用作对比
RS(k,j)=RS1(k,j)/2+RS2(k,j)/2;
span(j)=(length(xns_cal_peak)+length(xns_cal_peak1))/2*3e8/2/Laux2;
end

end


if(i==1)
    RS01(:,j)=RS1(:,j)
    RS02(:,j)=RS2(:,j)
    RS0(:,j)=RS(:,j)

end

% figure, hold on
% plot(xnfs)
% ylim([-2 2])
% ttmp=1:length(xnfs);
% plot(ttmp(mar_peak),xnfs(mar_peak),'ro')
% figure, hold on
% plot(xnfs1)
% ylim([-2 2])
% ttmp=1:length(xnfs1);
% plot(ttmp(mar_peak1),xnfs1(mar_peak1),'ro')

if(i==1)
sspe=20*log10(2*abs(fft(xns_cal_peak/length(xns_cal_peak))));
sspe=sspe(1:length(sspe)/2);
ffre=(0:length(sspe)-1)/length(sspe)*Laux1/2;
% figure,plot(ffre,sspe)
% set(0,'defaultfigurecolor','w')
% xlim([0 Laux1/2])
% set(gca,'FontSize',25,'FontName', 'Times New Roman')
% xlabel("Distance(m)","FontSize",25,'FontName', 'Times New Roman');
% ylabel('Power(dB)','FontSize',25,'FontName', 'Times New Roman')
else
sspe=20*log10(2*abs(fft(xns_cal_peak/length(xns_cal_peak))));
sspe=sspe(1:length(sspe)/2);
ffre=(0:length(sspe)-1)/length(sspe)*Laux2/2;
% figure,plot(ffre,sspe)
% set(0,'defaultfigurecolor','w')
% xlim([0 Laux2/2])
% set(gca,'FontSize',25,'FontName', 'Times New Roman')
% xlabel("Distance(m)","FontSize",25,'FontName', 'Times New Roman');
% ylabel('Power(dB)','FontSize',25,'FontName', 'Times New Roman')
end
% figure,plot(freqout,20*log10(2*abs(specout0/length(xns_cal_peak))))
% set(gca,'FontSize',25,'FontName', 'Times New Roman')
% xlabel("Distance(m)","FontSize",25,'FontName', 'Times New Roman');
% ylabel('Power(dB)','FontSize',25,'FontName', 'Times New Roman')

end
% figure,plot(freqout,20*log10(2*abs(specout0/length(xns_cal_peak))))
% set(gca,'FontSize',25,'FontName', 'Times New Roman')
% xlabel("Distance(m)","FontSize",25,'FontName', 'Times New Roman');
% ylabel('Power(dB)','FontSize',25,'FontName', 'Times New Roman')

%%

% std(((RS(:,j)*1000000)))
% % std(((RS1(:,j)*1000000)))
% % std(((RS2(:,j)*1000000)))
% mean(((RS(:,j))))
% std(((RS0(:,j)*1000000)))
% % std(((RS01(:,j)*1000000)))
% % std(((RS02(:,j)*1000000)))
% mean(((RS0(:,j))))

% figure,hold on
% plot(RS01(:,j))
% plot(RS02(:,j))
% plot(RS0(:,j))
% 
% figure,hold on
% plot(RS1(:,j))
% plot(RS2(:,j))
% plot(RS(:,j))

%%  


% fs=250e6;
% x=(1:length(xnfs))./fs;
% y=unwrap(angle(hilbert(xnfs)))./2./pi./(Laux2/3e8);
% 
% gamma1=max(y)/max(x);
% 
% para=polyfit(x,y,1);
% y0=para(1)*x+para(2);
% 
% RMSE1=(sum((y.'-y0).^2)./length(x)).^0.5/1e6;
% 
% figure,plot(x,y)
% 
% 
% x1=(1:length(xnfs1))./fs;
% y1=unwrap(angle(hilbert(xnfs1)))./2./pi./(Laux2/3e8);
% 
% gamma2=max(y1)/max(x1);
% para=polyfit(x1,y1,1);
% y01=para(1)*x1+para(2);
% 
% RMSE2=(sum((y1.'-y01).^2)./length(x1)).^0.5/1e6;
% 
%%
%%

    % --- 数据输入 ---
    x1 = Laux1;  y1 = mean(RS0(:,j));
    x2 = Laux2; y2 = mean(RS(:,j));
    if(abs(y1-y2)<0.1)
    res1=2*RS0;
    res2=2*RS0;
    else
    max_K = 150;    % 目标 K 值不超过 100
    tol = 0.5;      % K1 和 K2 允许的最大差值
    % ---------------
    
    fprintf('正在寻找 K <= %d，且 |K1 - K2| <= %.2f 的分支组合...\n', max_K, tol);
    fprintf('======================================================\n');
    
    % 根据 K=2mx±2y，估算 m 的最大搜索上限
    m1_max = ceil((max_K + 2*y1) / (2*x1));
    m2_max = ceil((max_K + 2*y2) / (2*x2));
    
    solutions_found = 0;
    index=[];
    
    for m1 = 0:m1_max
        for s1 = [-1, 1] % 点1：1为+，-1为-
            
            % 直接基于原始 x1, y1 计算 K1
            K1 = 2 * m1 * x1 + s1 * 2 * y1;
            
            % 过滤掉非物理域结果
            if K1 <= 0 || K1 > max_K
                continue;
            end
            
            for m2 = 1:m2_max
                for s2 = [-1, 1] % 点2：1为+，-1为-
                    
                    % 直接基于原始 x2, y2 计算 K2
                    K2 = 2 * m2 * x2 + s2 * 2 * y2;
                    
                    if K2 <= 0 || K2 > max_K
                        continue;
                    end
                    
                    % 核心判断：K1 和 K2 的差值是否在允许范围内
                    diff_K = abs(K1 - K2);
                    
                    if diff_K <= tol
                        solutions_found = solutions_found + 1;
                        
                        if s1 == 1, sign1 = '+'; else, sign1 = '-'; end
                        if s2 == 1, sign2 = '+'; else, sign2 = '-'; end
                        index=[index;m1,s1,m2,s2];
                        fprintf('找到可能组合 #%d (差值 %.4f):\n', solutions_found, diff_K);
                        fprintf('  点1: m1 = %d, 分支 %s\n', m1, sign1);
                        fprintf('       算得 K1 = 2*%d*x1 %s 2*y1 = %.6f\n', m1, sign1, K1);
                        fprintf('  点2: m2 = %d, 分支 %s\n', m2, sign2);
                        fprintf('       算得 K2 = 2*%d*x2 %s 2*y2 = %.6f\n', m2, sign2, K2);
                        fprintf('  > 推荐 K 值估计 (取均值) = %.6f\n', (K1 + K2)/2);
                        fprintf('------------------------------------------------------\n');
                    end
                end
            end
        end
    end
    
    if solutions_found == 0
        fprintf('未找到差值在 %.2f 以内的组合。\n', tol);
    else
        fprintf('共找到 %d 组符合条件的组合。\n', solutions_found);
    end

    group=1
    alpha=(index(group,4)*y2-index(group,2)*y1)./(index(group,1)*x1-index(group,3)*x2);

    res1=2*index(group,1)*alpha*x1+2*index(group,2)*RS0(:,j);
    res2=2*index(group,3)*alpha*x2+2*index(group,4)*RS(:,j);
    end
    figure,plot(res1/2)
    hold on
    plot(res2/2);
    set(0,'defaultfigurecolor','w')

    set(gca,'FontSize',25,'FontName', 'Times New Roman')
    xlabel("Index","FontSize",25,'FontName', 'Times New Roman');
    ylabel('Distance(m)','FontSize',25,'FontName', 'Times New Roman')
    res=[res1;res2];
%     figure,plot(res);
    final(j)=mean(res/2);
    var(j)=std(res/2)*1000000;
    final=final.';
    var=var.';
    std(res1/2)*1000000
    std(res2/2)*1000000
end
%%

figure,plot(final(tag,1),var(tag,1))
    
    set(0,'defaultfigurecolor','w')

    set(gca,'FontSize',25,'FontName', 'Times New Roman')
    xlabel("Distance(m)","FontSize",25,'FontName', 'Times New Roman');
    ylabel('Deviation(μm)','FontSize',25,'FontName', 'Times New Roman')