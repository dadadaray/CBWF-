%% 时间：2024/03/08
%% 功能：实现圆形vs直线定位+协同定位结果

clc;
clear all;
close all;
dbstop if error
%% 参数初始化
for repeat_time = 1:1 % 重复进行repeat_time次实验
    disp(repeat_time);
    sample_size = 1.6; % 定位区域离散化距离sample_size
    
    unit_length = 1.8; % 采样点1个单位间隔表示1.8米
    
    p_sim = 0.9; %相似度阈值
    
    randmon_num = 20; % 随机选择randmon_num个参考点用于计算圆心
    % 4个房间采集30秒的rp位置
    dir_train1 = 'train\';
    
    %% 区域1的初始及末尾顶点
    area1_ini = [0,0] * unit_length; area1_fin = [4.6,12] * unit_length; % 原来是4开始
    %% 区域2的初始及末尾顶点
    area2_ini = [3,5.89] * unit_length; area2_fin = [9.2,19.11] * unit_length; % 9.2 19.11
    %% 对上述两个区域进行边界扩张，写出定位具体的范围(5条障碍物)()
    position_boundary = [3.2,0,3.2,12.5; 3.2,5.5,0,5.5; 3.2,12.2,0,12.2;4.67,0,4.67,19.11;4.67,12.2,9.2,12.2;4.67,5.5,9.2,5.5] .* unit_length;
    
    pener_obstacle = 4; % 穿过一条障碍物补偿12dBm,因此-100dBm不准确，不带入计算，因此无需补偿
    AP_noise = 0;
    %% AP初始位置
    all_ap(1,:) = [2.7,9.33] + [randn,randn] * AP_noise;  all_ap(2,:) = [6.1,19.41] + [randn,randn] * AP_noise;  all_ap(3,:) = [1.6,6.67]+ [randn,randn] * AP_noise;
    
    all_ap(4,:) = [9.17,10.56]+ [randn,randn] * AP_noise; all_ap(5,:) = [6,12.22]+ [randn,randn] * AP_noise; all_ap(6,:) = [0.1,11.57]+ [randn,randn] * AP_noise;
    
    all_ap(7,:) = [2.8,1.83]+ [randn,randn] * AP_noise; all_ap(8,:) = [0.2,4.6]+ [randn,randn] * AP_noise; all_ap(9,:) = [1.67,-0.4]+ [randn,randn] * AP_noise;
    
    all_ap(10,:) = [9.1,19.11]+ [randn,randn] * AP_noise; all_ap(11,:) = [8.83,9.56]+ [randn,randn] * AP_noise;
    
    all_ap = all_ap.*unit_length; % AP的真实位置。
    
    %% 11个AP的mac地址
    all_ap_mac{1} = 'f0:b4:29:db:56:85'; all_ap_mac{2} = 'ec:60:73:e3:dd:14'; all_ap_mac{3} = 'ec:60:73:c7:50:8e';
    
    all_ap_mac{4} = 'f0:b4:29:da:ce:b6'; all_ap_mac{5} = '00:4b:f3:dc:59:46'; all_ap_mac{6} = '00:4b:f3:9f:ed:80';
    
    all_ap_mac{7} = 'b8:3a:08:b9:7f:d1'; all_ap_mac{8} = 'b8:3a:08:bb:43:11'; all_ap_mac{9} = 'b8:3a:08:ba:ca:31';
    
    all_ap_mac{10} = '00:4b:f3:ce:e8:1c'; all_ap_mac{11} = 'f0:b4:29:da:ce:6e';
    
    rand_rankap = randperm(size(all_ap,1)); % 产生1：size(all_ap,1)的随机数，且不重复
    
    for i = 1:11 % 所选择的AP数量
        ap_mac{i} = all_ap_mac{rand_rankap(i)};
        ap(i,:) = all_ap(rand_rankap(i),:);
    end
    
    %% 参考点的信息抽取
    [uni_allmac,mean_rssi,all_loc,mean_rssi_no_pen] = New_CreateFD_real(dir_train1,unit_length,ap_mac,pener_obstacle,position_boundary,ap,randmon_num);
    
    %% 找出所有可能配对的AP对，并记录他们，后续记录AP对之间的k值
    [same_x_y_ap] = Find_ap_coor(ap);
    
    %% 按照same_x_y_ap进行same_x_y_k值的确定
    [same_x_y_k,a_k,b_k,uncert_deg,Rsquare] = Cal_k_real(same_x_y_ap,mean_rssi,ap,all_loc);
    
    %% 根据same_x_y_k计算每一AP对所确定的圆circle_point_r,包含圆心及半径
    [circle_point_r] = Cal_cir_real(same_x_y_ap,same_x_y_k,ap);
    
    %% 为了方便计算，将整个定位区域进行离散化
    [vir_rp] = Gen_vir_real(area1_ini,area1_fin,area2_ini,area2_fin,sample_size);
    
    %% 确定每一个圆在离散点所划分的结果，按照same_x_y_k进行AP对的大小分配
    [vir_rp_ap_rank] = Gen_ap_rank_real(vir_rp,same_x_y_ap,circle_point_r,ap);
    
    %% 输入测试数据进行位置测试be
    dir_test = 'test\';
    %% 获得测试数据的原始信号强度
    [mean_rssi_test,all_loc_test] = New_CreateFD_test(dir_test,unit_length,ap_mac);
    
    %% 假设每一个测试点处于每一个虚拟点，则获得该点的信号补偿，进而获得该点的真实AP序列，最终得到定位结果
    tic;
    [test_loc,est_err,mis_rp_c,cc] = Test_est(vir_rp,ap,position_boundary,pener_obstacle,mean_rssi_test,vir_rp_ap_rank,all_loc_test,uncert_deg,p_sim,circle_point_r,Rsquare);
    t_cb = toc;
    
    %% 假设每一个测试点处于每一个虚拟点，则获得该点的信号补偿，进而获得该点的真实AP序列，最终得到定位结果(所有圆的权重均为1)
    %     [test_loc_no,est_err_no,mis_rp_c_no,cc_no] = Test_est_noweight(vir_rp,ap,position_boundary,pener_obstacle,mean_rssi_test,vir_rp_ap_rank,all_loc_test,uncert_deg,p_sim,circle_point_r,Rsquare);
    
    %%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% NN与KNN定位测试 %%%%%%%%%%%%%%%%%%%%%%%%%
    %% 进行定位
    tic
    [test_loc_NN,est_err_NN] = NN_estimate(mean_rssi_no_pen,mean_rssi_test,all_loc,all_loc_test);
    t_nn = toc;
    
    tic
    [test_loc_KNN,est_err_KNN] = KNN_estimate(mean_rssi_no_pen,mean_rssi_test,all_loc,all_loc_test);
    t_knn = toc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 直线边界定位 %%%%%%%%%%%%%%%%%%%%%%%%%
    %% 确定每一个圆在离散点所划分的结果，按照same_x_y_k进行AP对的大小分配
    [vir_rp_ap_line] = Gen_ap_rank_line(vir_rp,same_x_y_ap,ap);
    
    %% 假设每一个测试点处于每一个虚拟点，则获得该点的信号补偿，进而获得该点的真实AP序列，最终得到定位结果(直线)
    [test_loc_line,est_err_line,mis_rp_l,cc_line] = Test_est_line(vir_rp,ap,position_boundary,pener_obstacle,mean_rssi_test,vir_rp_ap_line,all_loc_test,p_sim,circle_point_r);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 下述为直线边界的协同定位结果进行定位 %%%%%%%%%%%%%%%%%%%%%%%%%
    %% 对于每一个测试点进行AP-rank，这里是利用原始信号强度进行排序，而非考虑障碍物
    [test_ap_rank] = Cal_test_rank(mean_rssi_test);
    
    %% 全部进行定位
    KK = 4; % 需要KK个点进行协同
    for W = 1:KK %选择多少个点来协同定位当前点
        for select_tp = 1:size(all_loc_test,1) % 当前定位select_tp点
            disp(select_tp);
            tp_index = [];
            %     [tp_index,top_similar] = Cal_ap_rank_sim(test_ap_rank,select_tp,mean_rssi_test,N,W);
            [tp_index,top_similar(W,select_tp)] = Cal_ap_rank_sim2(test_ap_rank,select_tp,mean_rssi_test,W);
            tp_index = [select_tp;tp_index]; % 把成对测试点的索引写进来
            collab_rssi = [];
            collab_rssi = mean_rssi_test(tp_index(1),:);
            for i = 2:length(tp_index)
                collab_rssi = [collab_rssi;mean_rssi_test(tp_index(i),:)];% 对应的成对测试点下接收信号强度
            end
            collab_cc = [];
            for i = 1:length(tp_index)
                collab_cc{i} = cc_line{tp_index(i)};% 这里确定是用cc还是cc_no，每一个定位点所可能位于的虚拟点
            end
            %% 基于所找到的最相似tp_index与此次测试点索引1进行peer-to-peer的协同定位
            tp_collab_true = all_loc_test(tp_index(1),:);
            for i = 2:length(tp_index)
                tp_collab_true = [tp_collab_true;all_loc_test(tp_index(i),:)]; % 真实位置
            end
            
            
            [collab_loc,collab_err] = Cal_collab(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true);  %对应的误差
            aa_line{W}(select_tp,:) = collab_loc(1,:); % 第W个协同点下的位置估计结果
            %     [collab_loc,collab_err] = Cal_collab_N(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true,N);  %N个点整体定位对应的误差
            %     [collab_loc,collab_err] = Cal_collab_logpath(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true);  %对应的误差
            co_err_line(W,select_tp) = collab_err(1); % 每一测试点的定位误差
        end
        mean_coerr_line(repeat_time,W) = mean(co_err_line(W,:)); %选择第W个相似点作为成对配对时的,定位误差
    end
    
    
    %% 位置协同结果
    fuse_loc_line = zeros(size(all_loc_test,1),2); % 定位点集结果初始化
    for i = 1:size(all_loc_test,1)
        for W = 1:KK
            fuse_loc_line(i,:) = fuse_loc_line(i,:) + aa_line{W}(i,:).*top_similar(W,i);
        end
        fuse_loc_line(i,:) = fuse_loc_line(i,:)/sum(top_similar(1:KK,i)); % 最终的估计结果
    end
    for i = 1:size(fuse_loc_line,1)
        err_fuse_line(i) = norm(fuse_loc_line(i,:)-all_loc_test(i,:));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 下述为圆形边界的协同定位结果进行定位 %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% 测试验证一下对于测试点select_tp的定位精度是否有提升
    %% 对于每一个测试点进行AP-rank，这里是利用原始信号强度进行排序，而非考虑障碍物
    tic;
    [test_ap_rank] = Cal_test_rank(mean_rssi_test);
    
    
    % cc = cc_no;
    %% 全部进行定位
    for W = 1:KK %选择多少个点来协同定位当前点
        for select_tp = 1:size(all_loc_test,1) % 当前定位select_tp点
            disp(select_tp);
            tp_index = [];
            %     [tp_index,top_similar] = Cal_ap_rank_sim(test_ap_rank,select_tp,mean_rssi_test,N,W);
            [tp_index,top_similar(W,select_tp)] = Cal_ap_rank_sim2(test_ap_rank,select_tp,mean_rssi_test,W);
            tp_index = [select_tp;tp_index]; % 把成对测试点的索引写进来
            collab_rssi = [];
            collab_rssi = mean_rssi_test(tp_index(1),:);
            for i = 2:length(tp_index)
                collab_rssi = [collab_rssi;mean_rssi_test(tp_index(i),:)];% 对应的成对测试点下接收信号强度
            end
            collab_cc = [];
            for i = 1:length(tp_index)
                collab_cc{i} = cc{tp_index(i)};% 这里确定是用cc还是cc_no，每一个定位点所可能位于的虚拟点
            end
            %% 基于所找到的最相似tp_index与此次测试点索引1进行peer-to-peer的协同定位
            tp_collab_true = all_loc_test(tp_index(1),:);
            for i = 2:length(tp_index)
                tp_collab_true = [tp_collab_true;all_loc_test(tp_index(i),:)]; % 真实位置
            end
            
            
            [collab_loc,collab_err] = Cal_collab(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true);  %对应的误差
            aa{W}(select_tp,:) = collab_loc(1,:); % 第W个协同点下的位置估计结果
            %     [collab_loc,collab_err] = Cal_collab_N(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true,N);  %对应的误差
            %     [collab_loc,collab_err] = Cal_collab_logpath(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true);  %对应的误差
            co_err(W,select_tp) = collab_err(1); % 每一测试点的定位误差
            
        end
        mean_coerr(repeat_time,W) = mean(co_err(W,:)); %选择第W个相似点作为成对配对时的,圆形边界定位误差
    end
    
    
    %% 位置协同结果
    fuse_loc = zeros(size(all_loc_test,1),2); % 定位点集结果初始化
    for i = 1:size(all_loc_test,1)
        for W = 1:KK
            fuse_loc(i,:) = fuse_loc(i,:) + aa{W}(i,:).*top_similar(W,i);
        end
        fuse_loc(i,:) = fuse_loc(i,:)/sum(top_similar(1:KK,i)); % 最终的估计结果
    end
    for i = 1:size(fuse_loc,1)
        err_fuse(i) = norm(fuse_loc(i,:)-all_loc_test(i,:));
    end
    t_cb2 = toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% OCLoc-N定位测试 %%%%%%%%%%%%%%%%%%%%%%%%%
%     tic
%     for select_tp = 1:size(all_loc_test,1)
%         disp(select_tp);
%         [oc_loc_N(select_tp,:),oc_err_N(select_tp)] = OCLoc_N(mean_rssi_no_pen,mean_rssi_test,all_loc,all_loc_test,KK,select_tp,test_ap_rank);
%     end
%     t_ocn = toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% OCLoc-N定位测试 %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% OCLoc-W定位测试 %%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    for select_tp = 1:size(all_loc_test,1)
        disp(select_tp);
        [oc_loc(select_tp,:),oc_err(select_tp)] = OCLoc(mean_rssi_no_pen,mean_rssi_test,all_loc,all_loc_test,KK,select_tp,test_ap_rank);
    end
    t_oc = toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% OCLoc-W定位测试 %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% 绘制CDF图
    %     figure;
    %     plot(sort(est_err_NN),[1/size(test_loc,1):1/size(test_loc,1):1],'k-.');
    %     lsr{1} = ['NN'];
    %     hold on
    %     plot(sort(est_err_KNN),[1/size(test_loc,1):1/size(test_loc,1):1],'g');
    %     lsr{2} = ['KNN'];
    %     hold on
    %     %%
    %     plot(sort(est_err_line),[1/125:1/125:1]);
    %     lsr{3} = ['Line'];
    %     hold on
    %     plot(sort(est_err),[1/125:1/125:1]);
    %     lsr{4} = ['CBWF'];
    %     hold on
    %     % for i = 1:KK
    %     %     plot(sort(co_err(i,:)),[1/125:1/125:1]);
    %     %     hold on
    %     %     lsr{i+4} = ['第',num2str(i)]; % 图注
    %     % end
    %     KK = 0;
    %     plot(sort(err_fuse_line),[1/125:1/125:1],':k');
    %     lsr{KK+5} = ['直线协同'];
    %     hold on
    %     plot(sort(err_fuse),[1/125:1/125:1],'r');
    %     lsr{KK+6} = ['圆形协同'];
    %     % legend('CBWF','最优','次优','季优','第四','多个');
    %     legend(lsr);
    %     xlabel('定位误差（米）');
    %     ylabel('CDF');
    %
    %     %% 输出结果
    %     fprintf('NN的误差为：%.2f\n',mean(est_err_NN));
    %     fprintf('KNN的误差为：%.2f\n',mean(est_err_KNN));
    %     fprintf('直线边界的误差为：%.2f米\n',mean(est_err_line));
    %     fprintf('圆形边界的误差为：%.2f米\n',mean(est_err));
    %     for i =1:KK
    %         fprintf('第%d个协同点下，优协同定位的误差为：%.2f米\n',i,mean_coerr(i));
    %     end
    %     fprintf('多个直线协同定位的误差为：%.2f米\n',mean(err_fuse_line));
    %     fprintf('多个圆形协同定位的误差为：%.2f米\n',mean(err_fuse));
    aa_NN(repeat_time) = mean(est_err_NN);
    aa_KNN(repeat_time) = mean(est_err_KNN);
    aa_li(repeat_time) = mean(est_err_line);
    aa_li_coll(repeat_time) = mean(err_fuse_line);
    aa_cb(repeat_time) = mean(est_err);
    aa_cb_coll(repeat_time) = mean(err_fuse);
    % 每一次循环下得到的平均定位结果
end
%% 绘制3D误差棒！
% figure;
% X = [3.87,3.51,3.75,3.12;3.60,3.11,3.36,2.85;3.53,3.09,3.03,2.62;3.65,3.34,3.02,2.71;3.86,3.78,3.51,3.47]; %(系列数据)
% E = [0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];%（标准差）
% legends = {'Line','Line\_col','CBWF','CBWF\_col'};%（图例）
% groupnames = cell(9,1);
% groupnames{1} = '1';groupnames{2} = '2';groupnames{3} = '3';groupnames{4} = '4';%（分组名称）
% groupnames{5} = '5';
% groupnames{6} = '5';groupnames{7} = '6';groupnames{8} = '7';%（分组名称）
% Title = '';
% Xlabel = '\tau';
% Ylabel = 'Localization accuracy (m)';
% barweb(X,E,1,groupnames,Title,Xlabel,Ylabel,jet,'none',legends,2,'plot');
% h= legend('Line','Line\_col','CBWF','CBWF\_col');
% box on;
% grid on;
% set(h,'Orientation','horizon');


%% 绘制CDF图
% figure;
% plot(sort(est_err_NN),[1/length(est_err_NN):1/length(est_err_NN):1],':r','LineWidth',1.5);
% hold on
% plot(sort(est_err_KNN),[1/length(est_err_NN):1/length(est_err_NN):1],'-.b','LineWidth',1.5);
% hold on
% plot(sort(est_err_line),[1/length(est_err_NN):1/length(est_err_NN):1],':k','LineWidth',1.5);
% hold on
% plot(sort(err_fuse_line),[1/length(est_err_NN):1/length(est_err_NN):1],'-.g','LineWidth',1.5);
% hold on
% plot(sort(est_err),[1/length(est_err_NN):1/length(est_err_NN):1],':m','LineWidth',1.5);
% hold on
% plot(sort(err_fuse),[1/length(est_err_NN):1/length(est_err_NN):1],'-.c','LineWidth',1.5);
% xlabel('Localization Accuracy(m)');
% ylabel('CDF');
% legend('NN','KNN','Line','Line\_col','CBWF','CBWF\_col','Location','SouthEast');
% 
% axis([0 8 0 1])



% close all;clear;clc;
% figure;
% X = [4.97,4.29,4.28;5.02,4.32,4.41;5.08,4.39,4.52;5.09,4.45,4.78;5.44,4.96,5.28;6.50,5.65,6.74;7.50,6.84,7.94;8.32,7.89,9.18]; %(系列数据)
% E = [0.09,0.07,0.13;0.12,0.12,0.21;0.24,0.18,0.18;0.22,0.23,0.61;0.58,0.38,0.50;0.86,0.49,1.98;1.29,1.25,1.88;1.65,2.03,2.71];%（标准差）
% groupnames = cell(8,1);
% groupnames{1} = '0.25';groupnames{2} = '0.5';groupnames{3} = '0.75';groupnames{4} = '1';%（分组名称）
% groupnames{5} = '2';groupnames{6} = '4';groupnames{7} = '6';groupnames{8} = '8';%（分组名称）
% Title = '';
% Xlabel = 'SD of Gaussian Noise';
% Ylabel = 'Localization accuracy (m)';
% barweb(X,E,1,groupnames,Title,Xlabel,Ylabel,jet,'none',legends,2,'plot');
% legend('Line','CBWF w/o \omega','CBWF w/ \omega');
% box on;
% grid on;