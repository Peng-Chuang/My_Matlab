% ---------------------------------------------------------------------------------
% Name	:  AFDD_Data_Convert_20190111.m  Helios02  单字节
% Ver	:  1.0
% Author:  MY
% Spec	:
% Env	:  Matlab2017a
% Idea	:
% His	:  ----------------------------------------------------------------------
%		   |	Date	|	Ver	  | 					Note
%   	   |------------|---------|----------------------------------------------
%		   |2019-01-11	|	1.0   | initial version
% ---------------------------------------------------------------------------------
% %
% % Existing Problems, list of not included functions:
% %  a)
% %	 b)

tic
clc;
clear all;



% 400M采样率
dir_name='E:\CNN\数据\400M采样率\L32A\PGA2\原始数据\临时存放2021';
% dir_name='E:\CNN\数据\400M采样率\L32A\PGA2\原始数据\x211207(方案三，去Y型采样，线圈两端采集dh前端并33nF电容串联2Ω电阻）';
% dir_name='E:\CNN\数据\400M采样率\L32A\PGA2\原始数据\x211209(方案三）';

% dir_name='E:\CNN\数据\400M采样率\L32A\PGA2\原始数据\220101(CW系列整机)';
% dir_name='E:\CNN\数据\400M采样率\L32A\PGA2\原始数据\220101申科模块方案';
% dir_name='E:\CNN\数据\400M采样率\L32A\PGA2\原始数据\直流电弧数据（申科互感器）';
dir_name='E:\CNN\数据\400M采样率\L32A\PGA2\原始数据\直流电弧数据（AFDD新整机）';
% dir_name='E:\CNN\数据\400M采样率\L32A\PGA2\原始数据\直流电弧自检数据（申科互感器）';
dir_name='E:\CNN\数据\400M采样率\L32A\PGA2\原始数据\220212(CW系列方案整机)';
% 50M采样率
% dir_name='E:\CNN\数据\50M采样率\L32A\PGA2\原始数据\x211210(方案三21dB）';
% dir_name='E:\CNN\数据\50M采样率\L32A\PGA2\原始数据\临时存放2021';
% dir_name='E:\CNN\数据\50M采样率\L32A\PGA2\原始数据\直流电弧数据（申科互感器）'
% dir_name='E:\CNN\数据\50M采样率\L32A\PGA2\原始数据\220101(CW系列整机)'
% dir_name='E:\CNN\数据\50M采样率\L32A\PGA2\原始数据\直流电弧数据（AFDD新整机）'

%% 通道个数
AVG_Ch_Num = 0;
ADC_Ch_Num = 0;

DISP_Ch_Num1 = 0;
FFT_Ch_Num1 = 0;
DISP_Ch_Num2 = 0;
FFT_Ch_Num2 = 0;

LF_I_Ch_Num = 0;
MSR_VI_Ch_Num = 0;
%% 通道选择数组
AVG_Ch_Sel = {};
ADC_Ch_Sel = {};

DISP_Ch_Sel1 = {};
FFT_Ch_Sel1 = {};

DISP_Ch_Sel2 = {};
FFT_Ch_Sel2 = {};

LF_I_Ch_Sel = {};
MSR_VI_Ch_Sel = {};
%% DISP channel Initial
DISP_Ch_Array1{1} = 'CH_AMP_DISP1';
DISP_Ch_Array1{2} = 'CH_T_DISP1';
DISP_Ch_Array1{3} = 'CH_PHA_JDG0';
DISP_Ch_Array1{4} = 'CH_AMP_JDG0';
DISP_Ch_Array1{5} = 'CHWAVE4M10M';
DISP_Ch_Array1{6} = 'CH_COR_JDG';

DISP_Ch_Array2{1} = 'CH_AMP_DISP0';
DISP_Ch_Array2{2} = 'CH_T_DISP0';
DISP_Ch_Array2{3} = 'CH_PHA_JDG1';
DISP_Ch_Array2{4} = 'CH_AMP_JDG1';
DISP_Ch_Array2{5} = 'CHWAVE8M48M';

% MSR_Ch_Array{1} = 'CH_MSR_V';
% MSR_Ch_Array{2} = 'CH_MSR_I';
% save -append 选项


if(regexp(dir_name,'400M'))



    fid=fopen('D:\工作学习资料\竞品测试\鼎信\第二版芯片\芯片采集程序\红外采集（马越）\20210521\20210520-10半波-400M_25_32A规格.ini');%数据截位1107-100半波-400M
    
%     fid=fopen('D:\工作学习资料\竞品测试\鼎信\第二版芯片\芯片采集程序\红外采集（马越）\20210521\20210520-2半波-400M_25_32A规格.ini');%数据截位1107-100半波-400M

    fid2=fopen('D:\工作学习资料\竞品测试\鼎信\第二版芯片\芯片采集程序\红外采集（马越）\20210521\20210520-100半波-400M_25_32A规格.ini');%数据截位1107-100半波-400M



    data = textscan(fid,'%s');
    IniParaCell1=data{1,1};

    data = textscan(fid2,'%s');
    IniParaCell2=data{1,1};

    IniParaCell=IniParaCell1;

else
    fid=fopen('D:\工作学习资料\竞品测试\鼎信\第二版芯片\芯片采集程序\红外采集（马越）\20210521\20210520-50半波-50M_红外.ini');%数据截位1107-100半波-400M
    
    fid2=fopen('D:\工作学习资料\竞品测试\鼎信\第二版芯片\芯片采集程序\红外采集（马越）\20210521\20210520-200半波-50M_网卡.ini');%数据截位1107-100半波-400M

    data = textscan(fid,'%s');
    IniParaCell1=data{1,1};

    data = textscan(fid2,'%s');
    IniParaCell2=data{1,1};

    IniParaCell=IniParaCell1;

end


if~isempty(regexp(IniParaCell{6}, '64', 'once'))%通道使能
    %% 通道特征值字节数


    AVG_Len = 64 + 4;   % PGA_Len = 4;
    ADC_Len = 51200;
    DISP_Len = 2;
    FFT_Len = 4;
    LF_I_Len = 4;
    MSR_VI_Len = 8;
elseif~isempty(regexp(IniParaCell{6}, '128', 'once'))
    AVG_Len = 64 + 4;   % PGA_Len = 4;
    ADC_Len = 51200;

    DISP_Len = 2;
    FFT_Len = 4;
    LF_I_Len = 4;
    MSR_VI_Len = 8;
end

for index = 35 : 1 :97
    if~isempty(regexp(IniParaCell{index}, 'true', 'once'))%通道使能
        if(index  == 35)  %幅值平均值
            AVG_Ch_Num = AVG_Ch_Num + 1;
            AVG_Ch_Sel{AVG_Ch_Num} = 'CHAVG';
        elseif(index <= 37)  %ADC
            ADC_Ch_Num = ADC_Ch_Num + 1;
            ADC_Ch_Sel{ADC_Ch_Num} = ['ADC',num2str(ADC_Ch_Num - 1)];
        elseif(index <= 43)  %高频离散度1
            DISP_Ch_Num1 = DISP_Ch_Num1 + 1;
            DISP_Ch_Sel1{DISP_Ch_Num1} = DISP_Ch_Array1{index - 37};

        elseif(index <= 52)  %FFT1
            FFT_Ch_Num1 = FFT_Ch_Num1 + 1;
            FFT_Ch_Sel1{FFT_Ch_Num1} = ['CH_FFT_',num2str(index-42),'M'];

        elseif(index <= 57)  %高频离散度2
            DISP_Ch_Num2 = DISP_Ch_Num2 + 1;
            DISP_Ch_Sel2{DISP_Ch_Num2} = DISP_Ch_Array2{index - 52};

        elseif(index <= 95)  %FFT2
            FFT_Ch_Num2 = FFT_Ch_Num2 + 1;
            FFT_Ch_Sel2{FFT_Ch_Num2} = ['CH_FFT_',num2str(index-47),'M'];

        elseif(index <= 96)  %低频电流
            continue;

        elseif(index <= 97)  %计量电压电流
            MSR_VI_Ch_Num = MSR_VI_Ch_Num + 1;
            continue;
        end

    end
end

%% 半波个数
BanBoNum = str2num(IniParaCell{5}(8:length(IniParaCell{5})));





%% 检查目录中是否存在已转换完成的matlab文件
fl=dirr(dir_name);
nl=length(fl);

cishu_mat=0;cishu_txt=0;
mat_name={};
txt_name={};
for i=1:nl
    chang(i)=length(fl(i).isdir);
    for j=1:chang(i)
        file = fullfile(fl(i).isdir(1).folder,fl(i).isdir(j).name);

        mingzi=fl(i).isdir(j).name;
        if(regexp(mingzi,'.mat'))
            cishu_mat=cishu_mat+1;
            mat_name{cishu_mat,1}=file(1:end-4);
            mat_name{cishu_mat,2}=mingzi(1:end-4);

        else
            cishu_txt=cishu_txt+1;
            txt_name{cishu_txt,1}=file(1:end-4);
            txt_name{cishu_txt,2}=mingzi(1:end-4);
        end
    end
end


if length(mat_name)==0

    file_name=txt_name(:,2);
    file_dir=txt_name(:,1);
else
    file_name=setdiff(txt_name(:,2),mat_name(:,2));
    file_dir=setdiff(txt_name(:,1),mat_name(:,1));
end



for i=1:length(file_name)
    cishu=i
    filename=[file_name{i},'.txt']


    if(regexp(filename,'网卡'))

        BanBoNum = str2num(IniParaCell2{5}(8:length(IniParaCell2{5})));
    end

    if(regexp(filename,'红外'))

        BanBoNum = str2num(IniParaCell1{5}(8:length(IniParaCell1{5})));
    end


    file=[file_dir{i},'.txt'];


    %     a=textread(file,'%s')';

    fid = fopen(file);
    data = textscan(fid,'%s');
    a=data{1,1};
    fclose(fid);


    %% 高频电弧特征值
    b = a;%高频电弧特征值
    Bytes_Length = length(b);
    valid_num=[];
    meter_valid_num=[];
    CHAVG=[];
    CH_PGA=[];
    ysAVG={};
    AVG_b = {};
    ADC_b = {};
    DISP1_b = {};
    FFT1_b = {};
    DISP2_b = {};
    FFT2_b = {};
    LF_I_b = {};
    MSR_VI_b = {};

    All_AVG_Bytes = 0;
    All_ADC_Bytes = 0;
    All_DISP_Bytes1 = 0;
    All_FFT_Bytes1 = 0;
    All_DISP_Bytes2 = 0;
    All_FFT_Bytes2 = 0;
    All_LF_I_Bytes = 0;
    All_MSR_VI_Bytes = 0;
    AllBytes = 0;
    %%  ----------------------------
    for k1 = 1 : BanBoNum
        %% 每个半波特征值个数列表
        b{end - 2 * BanBoNum + k1};
        valid_num(k1) = hex2dec(b{end - 2 * BanBoNum + k1});
        meter_valid_num(k1) = hex2dec(b{end - BanBoNum + k1});
        %% 不同通道特征值特征值个数
        AVG_Bytes	=   AVG_Len * AVG_Ch_Num;
        ADC_Bytes	=   ADC_Len * ADC_Ch_Num;

        DISP_Bytes1  = valid_num(k1) * DISP_Len * DISP_Ch_Num1;
        FFT_Bytes1   = valid_num(k1) * FFT_Len * FFT_Ch_Num1;
        DISP_Bytes2  = valid_num(k1) * DISP_Len * DISP_Ch_Num2;
        FFT_Bytes2   = valid_num(k1) * FFT_Len * FFT_Ch_Num2;
        LF_I_Bytes	=	valid_num(k1) * LF_I_Len * LF_I_Ch_Num;
        MSR_VI_Bytes	=	meter_valid_num(k1) * MSR_VI_Len * MSR_VI_Ch_Num;


        if(AVG_Bytes > 0)
            AVG_b(All_AVG_Bytes + 1 : All_AVG_Bytes + AVG_Bytes)  =  b(AllBytes + 1 : ...
                AllBytes + AVG_Bytes);
        end

        if(ADC_Bytes > 0)
            ADC_b(All_ADC_Bytes + 1 : All_ADC_Bytes + ADC_Bytes)  =  b(AllBytes + AVG_Bytes + 1 : ...
                AllBytes + AVG_Bytes + ADC_Bytes);
        end

        if(DISP_Bytes1 > 0)
            DISP1_b(All_DISP_Bytes1 + 1 : All_DISP_Bytes1 + DISP_Bytes1) =  b(AllBytes + AVG_Bytes + ADC_Bytes + 1 : ...
                AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1);
        end

        if(FFT_Bytes1 > 0)
            FFT1_b(All_FFT_Bytes1 + 1 : All_FFT_Bytes1 + FFT_Bytes1)  =  b(AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1 + 1 : ...
                AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1 + FFT_Bytes1);
        end


        if(DISP_Bytes2 > 0)
            DISP2_b(All_DISP_Bytes2 + 1 : All_DISP_Bytes2 + DISP_Bytes2) =  b(AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1 + FFT_Bytes1 + 1 : ...
                AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2);
        end


        if(FFT_Bytes2 > 0)
            FFT2_b(All_FFT_Bytes2 + 1 : All_FFT_Bytes2 + FFT_Bytes2) =  b(AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + 1 : ...
                AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2);
        end

        if(LF_I_Bytes > 0)
            LF_I_b(All_LF_I_Bytes + 1 : All_LF_I_Bytes + LF_I_Bytes) =  b(AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2 + 1 : ...
                AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2 + LF_I_Bytes);
        end

        if(MSR_VI_Bytes > 0)
            MSR_VI_b(All_MSR_VI_Bytes + 1 : All_MSR_VI_Bytes + MSR_VI_Bytes) =  b(AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2 + LF_I_Bytes + 1 : ...
                AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2 + LF_I_Bytes + MSR_VI_Bytes);
        end

        All_AVG_Bytes = All_AVG_Bytes + AVG_Bytes;
        All_ADC_Bytes = All_ADC_Bytes + ADC_Bytes;

        All_DISP_Bytes1 = All_DISP_Bytes1 + DISP_Bytes1;
        All_FFT_Bytes1 = All_FFT_Bytes1 + FFT_Bytes1;
        All_DISP_Bytes2 = All_DISP_Bytes2 + DISP_Bytes2;
        All_FFT_Bytes2 = All_FFT_Bytes2 + FFT_Bytes2;
        All_LF_I_Bytes = All_LF_I_Bytes + LF_I_Bytes;
        All_MSR_VI_Bytes = All_MSR_VI_Bytes + MSR_VI_Bytes;

        AllBytes = AllBytes + AVG_Bytes + ADC_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2 + LF_I_Bytes + MSR_VI_Bytes;
    end


    if(AVG_Bytes > 0)
        AVG_b_length = length(AVG_b);
        for i0  = 1: AVG_b_length/4
            ysAVG{i0} = [AVG_b{4*i0},AVG_b{4*i0 - 1},AVG_b{4*i0 - 2},AVG_b{4*i0 - 3}];

        end
        yAVG = hex2dec(ysAVG);
        yAVG=yAVG';
        for ii0 = 1:length(yAVG)/17
            CHAVG(16 * (ii0 - 1) + 1  :16 * (ii0 - 1) + 16) = yAVG(17 * (ii0 - 1) + 1  :17 * (ii0 - 1) + 16);
            CH_PGA(ii0) = yAVG(17 * (ii0 - 1) + 17) - 9;
        end

        if(exist([file(1:end-4),'.mat']))
            save([file(1:end-4),'.mat'],'CHAVG','-append');
        else
            save([file(1:end-4),'.mat'],'CHAVG');
        end
        %         end

        %% PGA
        if(exist([file(1:end-4),'.mat']))
            save([file(1:end-4),'.mat'],'CH_PGA','-append');
        else
            save([file(1:end-4),'.mat'],'CH_PGA');
        end
    end
    %             end
    %% ADC data   25600 51200bytes
    %------------------ADC数据解析待确认--------------------------------
    %     if(ADC_Bytes > 0)
    %         ADC_b_length = length(ADC_b);
    %         for i1  = 1: ADC_b_length/2
    %             ysADC{i1} = [ADC_b{2 * i1},ADC_b{2 * i1 - 1}];
    %             yADC(i0) = hex2dec(ysADC{i1});
    %         end
    %         for j = 1:length(yADC)
    %             if(yADC(j) > 2^11)
    %                 yADC(j) = mod(yADC(j),2^11) - 2^11;
    %             else
    %                 yADC(j) = yADC(j);
    %             end
    %         end
    %
    %         xADC = reshape(yADC,ADC_b_length/(2*BanBoNum),BanBoNum);
    %
    %         for l1 = 1:ADC_b_length/(ADC_Len*BanBoNum)
    %             for k1 = 1: BanBoNum
    %                 xADC1((k1-1)*(ADC_Len/2)+1:(k1-1)*(ADC_Len/2)+(ADC_Len/2),l1)  = xADC((l1-1)*(ADC_Len/2)+1:(l1-1)*(ADC_Len/2)+(ADC_Len/2),k1);
    %             end
    %         end
    %         xADC2 = xADC1';
    %
    %         ADC1 = xADC2(1,:);
    %         ADC2 = xADC2(2,:);
    %         if(exist([file(1:end-4),'.mat']))
    %             save([file(1:end-4),'.mat'],'ADC1','ADC2','-append');
    %         else
    %             save([file(1:end-4),'.mat'],'ADC1','ADC2');
    %         end
    %     end



    if(DISP_Bytes1 > 0)
        ysDISP1 = {};
        yDISP1 = [];
        DISP1_b_length = length(DISP1_b);


        DISP1_b_reshpae=reshape(DISP1_b,2,DISP1_b_length/2);
        yDISP1 = hex2dec(DISP1_b_reshpae);
        yDISP1_tmpppp=yDISP1(1:2:end)*1+yDISP1(2:2:end)*16^2;
        yDISP1=yDISP1_tmpppp';


        for ch = 1 : DISP_Ch_Num1
            eval([DISP_Ch_Sel1{ch},'=[];']);
        end
        AllDisp1ChNum = 0;
        for k1 = 1 : BanBoNum

            DISP1_s_Bytes_Tmp(k1) = valid_num(k1) * DISP_Len / 2;
            DISP1_Bytes_Tmp(k1)	=	DISP1_s_Bytes_Tmp(k1) * DISP_Ch_Num1 ;

            for ch = 1 : DISP_Ch_Num1
                eval([DISP_Ch_Sel1{ch},'= [',DISP_Ch_Sel1{ch},',  yDISP1(AllDisp1ChNum + DISP1_s_Bytes_Tmp(k1) * (ch - 1 ) + 1 : AllDisp1ChNum + DISP1_s_Bytes_Tmp(k1) * (ch - 1 ) + valid_num(k1))];']);
            end
            AllDisp1ChNum =  AllDisp1ChNum + DISP1_Bytes_Tmp(k1);
        end

        for ch = 1 : DISP_Ch_Num1
            if(exist([file(1:end-4),'.mat']))
                save([file(1:end-4),'.mat'],DISP_Ch_Sel1{ch},'-append');
            else
                save([file(1:end-4),'.mat'],DISP_Ch_Sel1{ch});
            end
        end

        if(exist('CH_COR_JDG'))
            cor_len = length(CH_COR_JDG);
            for CorNum = 1:cor_len
                if(CH_COR_JDG(CorNum) > 16)
                    CH_COR_JDG(CorNum) = -(CH_COR_JDG(CorNum) - 32 + 1);
                end
            end
            if(exist([file(1:end-4),'.mat']))
                save([file(1:end-4),'.mat'],'CH_COR_JDG','-append');
            else
                save([file(1:end-4),'.mat'],'CH_COR_JDG');
            end
        end


    end

    if(FFT_Bytes1 > 0)
        ysFFT1 = {};
        yFFT1 = [];
        FFT1_b_length = length(FFT1_b);

        FFT1_b_reshpae=reshape(FFT1_b,4,FFT1_b_length/4);
        yFFT1 = hex2dec(FFT1_b_reshpae);
        yFFT1_tmpppp=yFFT1(1:4:end)*1+yFFT1(2:4:end)*16^2+yFFT1(3:4:end)*16^4+yFFT1(4:4:end)*16^6;
        yFFT1=yFFT1_tmpppp';

        for ch = 1 : FFT_Ch_Num1
            eval([FFT_Ch_Sel1{ch},'=[];']);
        end
        AllFFT1ChNum = 0;
        for k1 = 1 : BanBoNum

            for ch = 1 : FFT_Ch_Num1
                eval([FFT_Ch_Sel1{ch},'= [',FFT_Ch_Sel1{ch},',  yFFT1(AllFFT1ChNum + valid_num(k1) * (ch - 1 ) + 1 : AllFFT1ChNum + valid_num(k1) * (ch - 1 ) + valid_num(k1))];']);
            end
            AllFFT1ChNum =  AllFFT1ChNum + valid_num(k1) * FFT_Ch_Num1;
        end


        for ch = 1 : FFT_Ch_Num1
            if(exist([file(1:end-4),'.mat']))
                save([file(1:end-4),'.mat'],FFT_Ch_Sel1{ch},'-append');
            else
                save([file(1:end-4),'.mat'],FFT_Ch_Sel1{ch});
            end
        end

    end

    if(DISP_Bytes2 > 0)
        ysDISP2 = {};
        yDISP2 = [];
        DISP2_b_length = length(DISP2_b);

        DISP2_b_reshpae=reshape(DISP2_b,2,DISP2_b_length/2);
        yDISP2 = hex2dec(DISP2_b_reshpae);
        yDISP2_tmpppp=yDISP2(1:2:end)*1+yDISP2(2:2:end)*16^2;
        yDISP2=yDISP2_tmpppp';


        for ch = 1 : DISP_Ch_Num2
            eval([DISP_Ch_Sel2{ch},'=[];']);
        end
        AllDisp2ChNum = 0;
        for k1 = 1 : BanBoNum

            DISP2_s_Bytes_Tmp(k1) = valid_num(k1) * DISP_Len / 2;
            DISP2_Bytes_Tmp(k1)	=	DISP2_s_Bytes_Tmp(k1) * DISP_Ch_Num2 ;

            for ch = 1 : DISP_Ch_Num2
                eval([DISP_Ch_Sel2{ch},'= [',DISP_Ch_Sel2{ch},',  yDISP2(AllDisp2ChNum + DISP2_s_Bytes_Tmp(k1) * (ch - 1 ) + 1 : AllDisp2ChNum + DISP2_s_Bytes_Tmp(k1) * (ch - 1 ) + valid_num(k1))];']);
            end
            AllDisp2ChNum =  AllDisp2ChNum + DISP2_Bytes_Tmp(k1);
        end

        for ch = 1 : DISP_Ch_Num2
            if(exist([file(1:end-4),'.mat']))
                save([file(1:end-4),'.mat'],DISP_Ch_Sel2{ch},'-append');
            else
                save([file(1:end-4),'.mat'],DISP_Ch_Sel2{ch});
            end
        end

    end


    if(FFT_Bytes2 > 0)
        ysFFT2 = {};
        yFFT2 = [];
        FFT2_b_length = length(FFT2_b);

        FFT2_b_reshpae=reshape(FFT2_b,4,FFT2_b_length/4);
        yFFT2 = hex2dec(FFT2_b_reshpae);
        yFFT2_tmpppp=yFFT2(1:4:end)*1+yFFT2(2:4:end)*16^2+yFFT2(3:4:end)*16^4+yFFT2(4:4:end)*16^6;
        yFFT2=yFFT2_tmpppp';


        for ch = 1 : FFT_Ch_Num2
            eval([FFT_Ch_Sel2{ch},'=[];']);
        end
        AllFFT2ChNum = 0;
        for k1 = 1 : BanBoNum

            for ch = 1 : FFT_Ch_Num2

                eval([FFT_Ch_Sel2{ch},'= [',FFT_Ch_Sel2{ch},',  yFFT2(AllFFT2ChNum + valid_num(k1) * (ch - 1 ) + 1 : AllFFT2ChNum + valid_num(k1) * (ch - 1 ) + valid_num(k1))];']);
            end
            AllFFT2ChNum =  AllFFT2ChNum + valid_num(k1) * FFT_Ch_Num2;
        end


        for ch = 1 : FFT_Ch_Num2
            if(exist([file(1:end-4),'.mat']))
                save([file(1:end-4),'.mat'],FFT_Ch_Sel2{ch},'-append');
            else
                save([file(1:end-4),'.mat'],FFT_Ch_Sel2{ch});
            end
        end

    end

    if(LF_I_Bytes > 0)
        ysLF_I = {};
        yLF_I = [];
        LF_I_b_length = length(LF_I_b);

        LF_I_b_reshpae=reshape(LF_I_b,4,LF_I_b_length/4);
        yLF_I = hex2dec(LF_I_b_reshpae);
        yLF_I_tmpppp=yLF_I(1:4:end)*1+yLF_I(2:4:end)*16^2+yLF_I(3:4:end)*16^4+yLF_I(4:4:end)*16^6;
        yLF_I=yLF_I_tmpppp';


        for j = 1:length(yLF_I)
            if(yLF_I(j) > 2^23)
                yLF_I(j) = mod(yLF_I(j),2^23) - 2^23;
            else
                yLF_I(j) = yLF_I(j);
            end
        end
        CHLF_I = yLF_I;
        if(exist([file(1:end-4),'.mat']))
            save([file(1:end-4),'.mat'],'CHLF_I','-append');
        else
            save([file(1:end-4),'.mat'],'CHLF_I');
        end
    end


    if(MSR_VI_Bytes > 0)
        ysMSR_VI = {};
        yMSR_VI = [];

        MSR_VI_b_length = length(MSR_VI_b);

        MSR_VI_b_reshpae=reshape(MSR_VI_b,4,MSR_VI_b_length/4);
        yMSR_VI = hex2dec(MSR_VI_b_reshpae);
        yMSR_VI_tmpppp=yMSR_VI(1:4:end)*1+yMSR_VI(2:4:end)*16^2+yMSR_VI(3:4:end)*16^4+yMSR_VI(4:4:end)*16^6;
        yMSR_VI=yMSR_VI_tmpppp';


        % 二进制补码装换
        for j = 1:length(yMSR_VI)
            if(yMSR_VI(j) > 2^23)
                yMSR_VI(j) = mod(yMSR_VI(j),2^23) - 2^23;
            else
                yMSR_VI(j) = yMSR_VI(j);
            end
        end
        %电压电流穿插解析
        CH_MSR_V =  yMSR_VI(1:2:end - 1);
        CH_MSR_I =  yMSR_VI(2:2:end);


        if(exist([file(1:end-4),'.mat']))
            save([file(1:end-4),'.mat'],'CH_MSR_V','CH_MSR_I','-append');
        else
            save([file(1:end-4),'.mat'],'CH_MSR_V','CH_MSR_I');
        end
    end

    save([file(1:end-4),'.mat'],'valid_num','meter_valid_num','-append');


end


function [list,sumbytes,varargout] = dirr(chemin,varargin)


verbose = 0;
% set to 1 to get folders list in command window

if nargin == 0
    chemin = cd;
end
if nargout == 0
    dum = varargin;
    varargin{1} = 'name';
    varargin = [varargin(1) dum];
end

fields = {'name' 'date' 'bytes' 'isdir'};

if regexp(chemin,'[\*\?]') % if chemin contains any ? or *
    filt = regexprep(chemin,'.*[\\/](.*\>)','$1');% get filter
    filt = regexprep(filt,'\.','\.');% in regexp format
    filt = regexprep(filt,'\*','.*');
    filt = regexprep(filt,'\?','.');
    filt = regexprep(filt,'(.*)','\\<$1');
    chemin = regexprep(chemin,'(.*)[\\/].*\>','$1');% and chemin
end

if not(isempty(varargin)) % if additional fields were provided after chemin
    for i = 1:length(fields)
        if strcmp(varargin{1},fields{i})% if first varargin matches a fieldname,
            % assume no filter was provided,

            if not(exist('filt','var'))% or it was in chemin and was set just before
                filt = '.*';% set it to wildcard
                break

            end
        end
    end
    if not(exist('filt','var'))% else
        filt = varargin{1};% first varargin is the filter
        varargin(1) = [];
    end
else% if no additional fields were provided and filter was not in chemin
    if not(exist('filt','var'))
        filt = '.*';
    end
end
% determine which varargin are fieldnames
whicharefields = zeros(1,length(varargin));
for i = 1:length(varargin)
    for j = 1:length(fields)
        if strcmp(varargin{i},fields{j})
            whicharefields(i) = 1;
            break
        end
    end
end
% set f2out and f2outfilt
f2out = {}; f2outfilt = {};
idx = 0;
if not(isempty(varargin))
    for i = 1:length(varargin)
        if whicharefields(i)
            idx = idx + 1;
            f2out{idx} = varargin{i};
            f2outfilt{idx} = '';
        else % if nargin{i} is not a fieldname, assume it's a filter
            f2outfilt{idx} = varargin{i};
        end
    end
end

%%%%%%%%%%%%%%%%%%%% START
if verbose
    disp(chemin);
end

list = dir(chemin);
if isempty(list)
    disp([chemin ' not found']);
    if nargout == 0
        clear list
    else
        for i = 1:nargout - 2
            varargout{i} = [];
        end
        sumbytes = 0;
    end
    return
end
% remove . and ..
i_file = 1;
while i_file <= length(list)
    if strcmp(list(i_file).name,'.')|strcmp(list(i_file).name,'..')
        list(i_file) = [];
    else
        i_file = i_file + 1;
    end
end

% set sumbytes
sumbytes = struct('total',0,'dir',{});
sumbytes(1).total = 0;
i_dir = 0;
% and all output fields
for i = 1:size(f2out,2)
    f2out{2,i} = {};
end
filenames = {};
todel = 0;
r = 1;
for i_out = 1:size(f2out,2)
    if strcmp(f2out{1,i_out},'isdir')
        if strcmp(f2outfilt{i_out},'0') % check if no recursion is wanted
            r = 0;
        end
    end
end

% for each item in list
for i_file = 1:length(list)
    for i_out = 1:size(f2out,2) % for every output field
        if not(isempty(f2outfilt{i_out}))% if there is a filter
            if strcmp(f2out{1,i_out},'bytes') % if field is 'bytes'
                line = [num2str(list(i_file).(f2out{1,i_out})) f2outfilt{i_out} ';']; % compare with filter numerically
                if eval(line)% if passes the filter
                    continue % continue to next field
                else
                    todel(end+1) = i_file; % else set to be deleted
                end
            elseif not(strcmp(f2out{1,i_out},'isdir'))% if field is 'name' or 'date'
                if regexpi(list(i_file).(f2out{1,i_out}),f2outfilt{i_out}) % apply filter
                    continue % continue to next field
                else
                    todel(end+1) = i_file; % else set to be deleted
                end
            end
        end
    end
    % once checked for every field's filter
    if todel(end) == i_file % if one didn't pass,
        if not(list(i_file).isdir) % and it's not a directory
            continue % skip this file and continue
        end
    else
        if regexpi(list(i_file).name,filt) % else, check for general filter on filename
            sumbytes(1).total = sumbytes(1).total + list(i_file).bytes; % sum bytes of that level
            for i_out = 1:size(f2out,2)% and assign all output fields with the values of that file
                f2out{2,i_out}{end+1} = [chemin filesep num2str(list(i_file).(f2out{1,i_out}))];
            end
        else
            todel(end+1) = i_file; % else the file will be removed from the list structure
        end
    end
    if list(i_file).isdir % if it's a directory
        if not(r)
            continue
        end
        i_dir = i_dir + 1;
        cheminext = strcat(chemin,filesep,list(i_file).name);
        % get it's content by recursion
        % write the line to enter eval
        line = '[list(i_file).isdir,sumbytes.dir(i_dir)';
        for i_out = 1:size(f2out,2)% with all the requested fields as temporary variables
            line = [line ',f2outtemp{' num2str(i_out) '}'];
        end
        line = [line '] = dirr(cheminext,filt'];
        for i_out = 1:size(f2out,2)
            line = [line ',f2out{1,' num2str(i_out) '}'];
            if f2outfilt{i_out}
                line = [line ',f2outfilt{' num2str(i_out) '}'];
            end
        end
        line = [line ');'];
        eval(line);

        for i_out = 1:size(f2out,2)
            f2out{2,i_out} = [f2out{2,i_out} f2outtemp{i_out}]; % catenate temporary variables with f2out
        end
        % sum bytes
        sumbytes(1).total = sumbytes(1).total + sumbytes(1).dir(i_dir).total; % that level + the next one
        list(i_file).bytes = sumbytes(1).dir(i_dir).total; % and set list(i_file).bytes to that value
        if list(i_file).bytes & todel(end) == i_file
            todel(end) = [];
        end
    end
end
todel(1) = [];
list(todel) = [];


for i_out = 1:size(f2out,2)
    varargout{i_out} = f2out{2,i_out};
end
if nargout == 0
    clear list
    disp(char(f2out{2,1}));
end
end
