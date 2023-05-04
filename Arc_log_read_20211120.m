% ---------------------------------------------------------------------------------
% Name	:  Arc_log_read.m 
% Ver	:  1.0
% Author:  JiSenrong
% Spec	:
% Env	:  Matlab2021a
% Idea	:
% His	:  ----------------------------------------------------------------------
%		   |	Date	|	Ver	  | 					Note
%   	   |------------|---------|----------------------------------------------
%		   |2021-11-08	|	1.0   | initial version
% ---------------------------------------------------------------------------------

clear;
clc;

% 通道特征值字节数
AVG_Len = 64;
DISP_Len = 2;
FFT_Len = 4;
MSR_V_Len = 4;
MSR_I_Len = 4;
  
% 通道个数
AVG_Ch_Num = 0;
DISP_Ch_Num1 = 0;
FFT_Ch_Num1 = 0;
DISP_Ch_Num2 = 0;
FFT_Ch_Num2 = 0;
MSR_V_Ch_Num = 0;
MSR_I_Ch_Num = 0;% 通道选择数组
AVG_Ch_Sel = {};
DISP_Ch_Sel1 = {};
FFT_Ch_Sel1 = {};

DISP_Ch_Sel2 = {};
FFT_Ch_Sel2 = {};

% DISP channel Initial
DISP_Ch_Array1{1} = 'CH_AMP_DISP1';
DISP_Ch_Array1{2} = 'CH_T_DISP1';
DISP_Ch_Array1{3} = 'CH_PHA_JDG0';
DISP_Ch_Array1{4} = 'CH_AMP_JDG0';
DISP_Ch_Array1{5} = 'CHWAVE4M10M';
DISP_Ch_Array1{6} = 'CH_COR_JDG';

DISP_Ch_Array2{1} = 'CH_AMP_DISP0';
DISP_Ch_Array2{2} = 'CH_T_DISP0';
DISP_Ch_Array2{3} = 'CH_PHA_JDG1';
DISP_Ch_Array2{4} = 'CHWAVE8M48M';

CH_PGA=-9;

% 各通道名称
for index = 1:60
    if(index  == 1)  %幅值平均值
        AVG_Ch_Num = AVG_Ch_Num + 1;
        AVG_Ch_Sel{AVG_Ch_Num} = 'CHAVG';
    elseif(index <= 7)   %高频离散度1
        DISP_Ch_Num1 = DISP_Ch_Num1 + 1;
        DISP_Ch_Sel1{DISP_Ch_Num1} = DISP_Ch_Array1{index-1};
    elseif(index <= 16)  % FFT1
        FFT_Ch_Num1 = FFT_Ch_Num1 + 1;
        FFT_Ch_Sel1{FFT_Ch_Num1} = ['CH_FFT_',num2str(index-6),'M'];
    elseif(index <= 20)  %高频离散度2
        DISP_Ch_Num2 = DISP_Ch_Num2 + 1;
        DISP_Ch_Sel2{DISP_Ch_Num2} = DISP_Ch_Array2{index - 16};
    elseif(index <= 58)  %FFT2
        FFT_Ch_Num2 = FFT_Ch_Num2 + 1;
        FFT_Ch_Sel2{FFT_Ch_Num2} = ['CH_FFT_',num2str(index-10),'M'];
    elseif(index <= 59)  %计量电压
        MSR_V_Ch_Num = MSR_V_Ch_Num + 1;
        continue;
    else    %计量电流
        MSR_I_Ch_Num = MSR_I_Ch_Num + 1;
        continue;
    end
end

% -------------------------------------------------- 文件路径 ---------------------------------------- %
dir_name = 'D:\Topscomm\AFDD电弧日志解析\第一版\0615';
% dir_name = 'D:\Topscomm\AFDD软件调试助手\Afdd_Assistant_V22.05.21\Afdd_Assistant_V22.04.27';

fl=dirr(dir_name);
nl=length(fl);

cishu_mat=0;
cishu_txt=0;
mat_name={};
txt_name={};

for i=1:nl
    chang(i)=length(fl(i).isdir);
    for j=1:chang(i)
        
        %%%%%%%%%%%%%%%%%%% 适用于多个子文件夹 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        file = fullfile(fl(i).isdir(1).folder,fl(i).isdir(j).name);
        mingzi=fl(i).isdir(j).name;

        
        %%%%%%%%%%%%%%%%%%% 适用于无子文件夹 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         file = fullfile(fl(i).folder,fl(i).name);
%         mingzi=fl(i).name;
        
        
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
    
    file=[file_dir{i},'.txt'];
    fid = fopen(file, 'r');
    ParaCell = textscan(fid, '%s');
    all_content = ParaCell{1,1};
    fclose(fid);
    
    % --------------------------------------- 删除除电弧特征外的多余数据 ----------------------------------- %
    index_flag = find(strcmp('设备描述信息：智能开关研发部-智能开关', all_content));
    
    % 如果不清空，zhwen_bytes_index/valid_bytes_index 会导致本次错误。 
    % 例如，前一次4个半波，本次1个半波，则上次的第1个半波被覆盖，但上次的后3个半波仍存有值。
    zhwen_bytes_index = [];            
    valid_bytes_index = [];
    
    for k=1:length(index_flag)
        valid_num(k) = hex2dec(all_content{index_flag(k)+1});
        meter_valid_num(k) = hex2dec(all_content{index_flag(k)+5});
        zhwen_bytes_index(6*k-5:6*k,1) = (index_flag(k)-5:index_flag(k))';
        valid_bytes_index(8*k-7:8*k,1) = (index_flag(k)+1:index_flag(k)+8)';
    end
    
    begin_del_num = 18 * 4;
    del_patch_index = index_flag + begin_del_num;
    del_index_patch = [];
    
    for ii=1:length(del_patch_index)
        for jj=1:59
            if (jj <= 6)
                base = valid_num(ii) * DISP_Len;
                data_end = ceil(base/20) * 20;
                differ = data_end - base;
                del_index = del_patch_index(ii) + data_end * jj -differ+1:del_patch_index(ii) + data_end * jj;
            elseif (jj<=15)
                base = valid_num(ii) * FFT_Len;
                data_end = ceil(base/20) * 20;
                differ = data_end - base;
                del_index = del_patch_index(ii) + data_end * (jj-3) -differ+1:del_patch_index(ii) + data_end * (jj-3); % jj-6+(6/2)=jj-3
            elseif (jj<=19)
                base = valid_num(ii) * DISP_Len;
                data_end = ceil(base/20) * 20;
                differ_1 = data_end - base;
                del_index = del_patch_index(ii) + data_end * (jj+9) -differ_1+1:del_patch_index(ii) + data_end * (jj+9); % jj-15+6+2*9=jj+9
            elseif (jj<=57)
                base = valid_num(ii) * FFT_Len;
                data_end = ceil(base/20) * 20;
                differ_2 = data_end - base;
                del_index = del_patch_index(ii) + data_end * (jj-5) -differ_2+1:del_patch_index(ii) + data_end * (jj-5); % jj-19+(6/2)+9+(4/2)=jj-5
            elseif (jj<=59)
                base = meter_valid_num(ii) * MSR_V_Len;
                data_end = ceil(base/20) * 20;
                differ = data_end - base;
                if differ_2==0  % 特殊情况，47个64位通道不需要补0
                   del_index = del_patch_index(ii) + data_end * (jj-5)-47*20-differ+1:del_patch_index(ii) + data_end * (jj-5)-47*20;
%                 elseif differ_1==0  % 特殊情况，10个32位通道的FFT不需要补0
%                     ...                   
                else
                    del_index = del_patch_index(ii) + data_end * (jj-5) -differ+1:del_patch_index(ii) + data_end * (jj-5); % jj-57+(6/2)+9+(4/2)+38=jj-5 
                end
            end
            del_index_patch = [del_index_patch, del_index];
        end
    end
    
    del_index_all = [zhwen_bytes_index;valid_bytes_index;del_index_patch']; % 需删除的所有索引
    refer = all_content;
    
    all_content(del_index_all) = [];
    
    % 删除0x
    for kk=1:length(all_content)
        temp=all_content{kk};
        data_dh{kk,1}=temp(3:4);
    end
    
    % ------------------------------------------------ 解析16进制数据 ----------------------------------------------- %
    % 半波个数
    BanBoNum = length(valid_num);
    
    % 高频电弧特征值
    b = data_dh;
    Bytes_Length = length(b);
    AVG_b = {};
    ADC_b = {};
    DISP1_b = {};
    FFT1_b = {};
    DISP2_b = {};
    FFT2_b = {};
    MSR_V_b = {};
    MSR_I_b = {};
    
    All_AVG_Bytes = 0;
    All_DISP_Bytes1 = 0;
    All_FFT_Bytes1 = 0;
    All_DISP_Bytes2 = 0;
    All_FFT_Bytes2 = 0;
    All_MSR_V_Bytes = 0;
    All_MSR_I_Bytes = 0;
    AllBytes = 0;
    
    for k1 = 1 : BanBoNum
        % 不同通道特征值特征值个数
        AVG_Bytes	=   AVG_Len * AVG_Ch_Num;
        DISP_Bytes1  = valid_num(k1) * DISP_Len * DISP_Ch_Num1;
        FFT_Bytes1   = valid_num(k1) * FFT_Len * FFT_Ch_Num1;
        DISP_Bytes2  = valid_num(k1) * DISP_Len * DISP_Ch_Num2;
        FFT_Bytes2   = valid_num(k1) * FFT_Len * FFT_Ch_Num2;
        MSR_V_Bytes = meter_valid_num(k1) * MSR_V_Len * MSR_V_Ch_Num;
        MSR_I_Bytes = meter_valid_num(k1) * MSR_I_Len * MSR_I_Ch_Num;
        
        if(AVG_Bytes > 0)
            AVG_b(All_AVG_Bytes + 1 : All_AVG_Bytes + AVG_Bytes)  =  b(AllBytes + 1 : ...
                AllBytes + AVG_Bytes);
        end
        
        if(DISP_Bytes1 > 0)
            DISP1_b(All_DISP_Bytes1 + 1 : All_DISP_Bytes1 + DISP_Bytes1) =  b(AllBytes + AVG_Bytes + 1: AllBytes + AVG_Bytes + DISP_Bytes1);
        end
        
        if(FFT_Bytes1 > 0)
            FFT1_b(All_FFT_Bytes1 + 1 : All_FFT_Bytes1 + FFT_Bytes1)  =  b(AllBytes + AVG_Bytes + DISP_Bytes1 + 1 : ...
                AllBytes + AVG_Bytes + DISP_Bytes1 + FFT_Bytes1);
        end
        
        if(DISP_Bytes2 > 0)
            DISP2_b(All_DISP_Bytes2 + 1 : All_DISP_Bytes2 + DISP_Bytes2) =  b(AllBytes + AVG_Bytes + DISP_Bytes1 + FFT_Bytes1 + 1 : ...
                AllBytes + AVG_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2);
        end
        
        if(FFT_Bytes2 > 0)
            FFT2_b(All_FFT_Bytes2 + 1 : All_FFT_Bytes2 + FFT_Bytes2) =  b(AllBytes + AVG_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + 1 : ...
                AllBytes + AVG_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2);
        end
        
        if(MSR_V_Bytes > 0)
            MSR_V_b(All_MSR_V_Bytes + 1 : All_MSR_V_Bytes + MSR_V_Bytes) =  b(AllBytes + AVG_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2 + 1 : ...
                AllBytes + AVG_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2 + MSR_V_Bytes);
        end
        
        if(MSR_I_Bytes > 0)
            MSR_I_b(All_MSR_I_Bytes + 1 : All_MSR_I_Bytes + MSR_I_Bytes) =  b(AllBytes + AVG_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2 + MSR_V_Bytes + 1 : ...
                AllBytes + AVG_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2 + MSR_V_Bytes + MSR_I_Bytes);
        end
        
        All_AVG_Bytes = All_AVG_Bytes + AVG_Bytes;
        All_DISP_Bytes1 = All_DISP_Bytes1 + DISP_Bytes1;
        All_FFT_Bytes1 = All_FFT_Bytes1 + FFT_Bytes1;
        All_DISP_Bytes2 = All_DISP_Bytes2 + DISP_Bytes2;
        All_FFT_Bytes2 = All_FFT_Bytes2 + FFT_Bytes2;
        All_MSR_V_Bytes = All_MSR_V_Bytes + MSR_V_Bytes;
        All_MSR_I_Bytes = All_MSR_I_Bytes + MSR_I_Bytes;
        
        AllBytes = AllBytes + AVG_Bytes + DISP_Bytes1 + FFT_Bytes1 + DISP_Bytes2 + FFT_Bytes2 + MSR_V_Bytes + MSR_I_Bytes;
    end
    
    if(AVG_Bytes > 0)
        AVG_b_length = length(AVG_b);
        for i0  = 1: AVG_b_length/4
            ysAVG{i0} = [AVG_b{4*i0},AVG_b{4*i0 - 1},AVG_b{4*i0 - 2},AVG_b{4*i0 - 3}];
        end
        yAVG = hex2dec(ysAVG);
        CHAVG=yAVG';
        
        if(exist([file(1:end-4),'.mat']))
            save([file(1:end-4),'.mat'],'CHAVG','-append');
        else
            save([file(1:end-4),'.mat'],'CHAVG');
        end
    end
    
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
                eval([FFT_Ch_Sel1{ch},'= [',FFT_Ch_Sel1{ch},', yFFT1(AllFFT1ChNum + valid_num(k1) * (ch - 1 ) + 1 : AllFFT1ChNum + valid_num(k1) * (ch - 1 ) + valid_num(k1))];']);
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
    
    if(MSR_V_Bytes > 0)
        ysMSR_V = {};
        yMSR_V = [];
        
        MSR_V_b_length = length(MSR_V_b);
        
        MSR_V_b_reshpae=reshape(MSR_V_b,4,MSR_V_b_length/4);
        yMSR_V = hex2dec(MSR_V_b_reshpae);
        yMSR_V_tmpppp=yMSR_V(1:4:end)*1+yMSR_V(2:4:end)*16^2+yMSR_V(3:4:end)*16^4+yMSR_V(4:4:end)*16^6;
        yMSR_V=yMSR_V_tmpppp';
        
        % 二进制补码装换
        for j = 1:length(yMSR_V)
            if(yMSR_V(j) > 2^23)
                yMSR_V(j) = mod(yMSR_V(j),2^23) - 2^23;
            else
                yMSR_V(j) = yMSR_V(j);
            end
        end
        %电压解析
        CH_MSR_V =  yMSR_V;
        
        if(exist([file(1:end-4),'.mat']))
            save([file(1:end-4),'.mat'],'CH_MSR_V','-append');
        else
            save([file(1:end-4),'.mat'],'CH_MSR_V');
        end
    end
    
    if(MSR_I_Bytes > 0)
        ysMSR_I = {};
        yMSR_I = [];
        
        MSR_I_b_length = length(MSR_I_b);
        
        MSR_I_b_reshpae=reshape(MSR_I_b,4,MSR_I_b_length/4);
        yMSR_I = hex2dec(MSR_I_b_reshpae);
        yMSR_I_tmpppp=yMSR_I(1:4:end)*1+yMSR_I(2:4:end)*16^2+yMSR_I(3:4:end)*16^4+yMSR_I(4:4:end)*16^6;
        yMSR_I=yMSR_I_tmpppp';
        
        % 二进制补码装换
        for j = 1:length(yMSR_I)
            if(yMSR_I(j) > 2^23)
                yMSR_I(j) = mod(yMSR_I(j),2^23) - 2^23;
            else
                yMSR_I(j) = yMSR_I(j);
            end
        end
        % 电流解析
        CH_MSR_I =  yMSR_I;
        
        if(exist([file(1:end-4),'.mat']))
            save([file(1:end-4),'.mat'],'CH_MSR_I','-append');
        else
            save([file(1:end-4),'.mat'],'CH_MSR_I');
        end
    end
    
    save([file(1:end-4),'.mat'],'valid_num','meter_valid_num','CH_PGA','-append');
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


