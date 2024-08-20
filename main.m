%% ʱ�䣺2024/03/08
%% ���ܣ�ʵ��Բ��vsֱ�߶�λ+Эͬ��λ���

clc;
clear all;
close all;
dbstop if error
%% ������ʼ��
for repeat_time = 1:1 % �ظ�����repeat_time��ʵ��
    disp(repeat_time);
    sample_size = 1.6; % ��λ������ɢ������sample_size
    
    unit_length = 1.8; % ������1����λ�����ʾ1.8��
    
    p_sim = 0.9; %���ƶ���ֵ
    
    randmon_num = 20; % ���ѡ��randmon_num���ο������ڼ���Բ��
    % 4������ɼ�30���rpλ��
    dir_train1 = 'train\';
    
    %% ����1�ĳ�ʼ��ĩβ����
    area1_ini = [0,0] * unit_length; area1_fin = [4.6,12] * unit_length; % ԭ����4��ʼ
    %% ����2�ĳ�ʼ��ĩβ����
    area2_ini = [3,5.89] * unit_length; area2_fin = [9.2,19.11] * unit_length; % 9.2 19.11
    %% ����������������б߽����ţ�д����λ����ķ�Χ(5���ϰ���)()
    position_boundary = [3.2,0,3.2,12.5; 3.2,5.5,0,5.5; 3.2,12.2,0,12.2;4.67,0,4.67,19.11;4.67,12.2,9.2,12.2;4.67,5.5,9.2,5.5] .* unit_length;
    
    pener_obstacle = 4; % ����һ���ϰ��ﲹ��12dBm,���-100dBm��׼ȷ����������㣬������貹��
    AP_noise = 0;
    %% AP��ʼλ��
    all_ap(1,:) = [2.7,9.33] + [randn,randn] * AP_noise;  all_ap(2,:) = [6.1,19.41] + [randn,randn] * AP_noise;  all_ap(3,:) = [1.6,6.67]+ [randn,randn] * AP_noise;
    
    all_ap(4,:) = [9.17,10.56]+ [randn,randn] * AP_noise; all_ap(5,:) = [6,12.22]+ [randn,randn] * AP_noise; all_ap(6,:) = [0.1,11.57]+ [randn,randn] * AP_noise;
    
    all_ap(7,:) = [2.8,1.83]+ [randn,randn] * AP_noise; all_ap(8,:) = [0.2,4.6]+ [randn,randn] * AP_noise; all_ap(9,:) = [1.67,-0.4]+ [randn,randn] * AP_noise;
    
    all_ap(10,:) = [9.1,19.11]+ [randn,randn] * AP_noise; all_ap(11,:) = [8.83,9.56]+ [randn,randn] * AP_noise;
    
    all_ap = all_ap.*unit_length; % AP����ʵλ�á�
    
    %% 11��AP��mac��ַ
    all_ap_mac{1} = 'f0:b4:29:db:56:85'; all_ap_mac{2} = 'ec:60:73:e3:dd:14'; all_ap_mac{3} = 'ec:60:73:c7:50:8e';
    
    all_ap_mac{4} = 'f0:b4:29:da:ce:b6'; all_ap_mac{5} = '00:4b:f3:dc:59:46'; all_ap_mac{6} = '00:4b:f3:9f:ed:80';
    
    all_ap_mac{7} = 'b8:3a:08:b9:7f:d1'; all_ap_mac{8} = 'b8:3a:08:bb:43:11'; all_ap_mac{9} = 'b8:3a:08:ba:ca:31';
    
    all_ap_mac{10} = '00:4b:f3:ce:e8:1c'; all_ap_mac{11} = 'f0:b4:29:da:ce:6e';
    
    rand_rankap = randperm(size(all_ap,1)); % ����1��size(all_ap,1)����������Ҳ��ظ�
    
    for i = 1:11 % ��ѡ���AP����
        ap_mac{i} = all_ap_mac{rand_rankap(i)};
        ap(i,:) = all_ap(rand_rankap(i),:);
    end
    
    %% �ο������Ϣ��ȡ
    [uni_allmac,mean_rssi,all_loc,mean_rssi_no_pen] = New_CreateFD_real(dir_train1,unit_length,ap_mac,pener_obstacle,position_boundary,ap,randmon_num);
    
    %% �ҳ����п�����Ե�AP�ԣ�����¼���ǣ�������¼AP��֮���kֵ
    [same_x_y_ap] = Find_ap_coor(ap);
    
    %% ����same_x_y_ap����same_x_y_kֵ��ȷ��
    [same_x_y_k,a_k,b_k,uncert_deg,Rsquare] = Cal_k_real(same_x_y_ap,mean_rssi,ap,all_loc);
    
    %% ����same_x_y_k����ÿһAP����ȷ����Բcircle_point_r,����Բ�ļ��뾶
    [circle_point_r] = Cal_cir_real(same_x_y_ap,same_x_y_k,ap);
    
    %% Ϊ�˷�����㣬��������λ���������ɢ��
    [vir_rp] = Gen_vir_real(area1_ini,area1_fin,area2_ini,area2_fin,sample_size);
    
    %% ȷ��ÿһ��Բ����ɢ�������ֵĽ��������same_x_y_k����AP�ԵĴ�С����
    [vir_rp_ap_rank] = Gen_ap_rank_real(vir_rp,same_x_y_ap,circle_point_r,ap);
    
    %% ����������ݽ���λ�ò���be
    dir_test = 'test\';
    %% ��ò������ݵ�ԭʼ�ź�ǿ��
    [mean_rssi_test,all_loc_test] = New_CreateFD_test(dir_test,unit_length,ap_mac);
    
    %% ����ÿһ�����Ե㴦��ÿһ������㣬���øõ���źŲ�����������øõ����ʵAP���У����յõ���λ���
    tic;
    [test_loc,est_err,mis_rp_c,cc] = Test_est(vir_rp,ap,position_boundary,pener_obstacle,mean_rssi_test,vir_rp_ap_rank,all_loc_test,uncert_deg,p_sim,circle_point_r,Rsquare);
    t_cb = toc;
    
    %% ����ÿһ�����Ե㴦��ÿһ������㣬���øõ���źŲ�����������øõ����ʵAP���У����յõ���λ���(����Բ��Ȩ�ؾ�Ϊ1)
    %     [test_loc_no,est_err_no,mis_rp_c_no,cc_no] = Test_est_noweight(vir_rp,ap,position_boundary,pener_obstacle,mean_rssi_test,vir_rp_ap_rank,all_loc_test,uncert_deg,p_sim,circle_point_r,Rsquare);
    
    %%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% NN��KNN��λ���� %%%%%%%%%%%%%%%%%%%%%%%%%
    %% ���ж�λ
    tic
    [test_loc_NN,est_err_NN] = NN_estimate(mean_rssi_no_pen,mean_rssi_test,all_loc,all_loc_test);
    t_nn = toc;
    
    tic
    [test_loc_KNN,est_err_KNN] = KNN_estimate(mean_rssi_no_pen,mean_rssi_test,all_loc,all_loc_test);
    t_knn = toc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% ֱ�߽߱綨λ %%%%%%%%%%%%%%%%%%%%%%%%%
    %% ȷ��ÿһ��Բ����ɢ�������ֵĽ��������same_x_y_k����AP�ԵĴ�С����
    [vir_rp_ap_line] = Gen_ap_rank_line(vir_rp,same_x_y_ap,ap);
    
    %% ����ÿһ�����Ե㴦��ÿһ������㣬���øõ���źŲ�����������øõ����ʵAP���У����յõ���λ���(ֱ��)
    [test_loc_line,est_err_line,mis_rp_l,cc_line] = Test_est_line(vir_rp,ap,position_boundary,pener_obstacle,mean_rssi_test,vir_rp_ap_line,all_loc_test,p_sim,circle_point_r);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% ����Ϊֱ�߽߱��Эͬ��λ������ж�λ %%%%%%%%%%%%%%%%%%%%%%%%%
    %% ����ÿһ�����Ե����AP-rank������������ԭʼ�ź�ǿ�Ƚ������򣬶��ǿ����ϰ���
    [test_ap_rank] = Cal_test_rank(mean_rssi_test);
    
    %% ȫ�����ж�λ
    KK = 4; % ��ҪKK�������Эͬ
    for W = 1:KK %ѡ����ٸ�����Эͬ��λ��ǰ��
        for select_tp = 1:size(all_loc_test,1) % ��ǰ��λselect_tp��
            disp(select_tp);
            tp_index = [];
            %     [tp_index,top_similar] = Cal_ap_rank_sim(test_ap_rank,select_tp,mean_rssi_test,N,W);
            [tp_index,top_similar(W,select_tp)] = Cal_ap_rank_sim2(test_ap_rank,select_tp,mean_rssi_test,W);
            tp_index = [select_tp;tp_index]; % �ѳɶԲ��Ե������д����
            collab_rssi = [];
            collab_rssi = mean_rssi_test(tp_index(1),:);
            for i = 2:length(tp_index)
                collab_rssi = [collab_rssi;mean_rssi_test(tp_index(i),:)];% ��Ӧ�ĳɶԲ��Ե��½����ź�ǿ��
            end
            collab_cc = [];
            for i = 1:length(tp_index)
                collab_cc{i} = cc_line{tp_index(i)};% ����ȷ������cc����cc_no��ÿһ����λ��������λ�ڵ������
            end
            %% �������ҵ���������tp_index��˴β��Ե�����1����peer-to-peer��Эͬ��λ
            tp_collab_true = all_loc_test(tp_index(1),:);
            for i = 2:length(tp_index)
                tp_collab_true = [tp_collab_true;all_loc_test(tp_index(i),:)]; % ��ʵλ��
            end
            
            
            [collab_loc,collab_err] = Cal_collab(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true);  %��Ӧ�����
            aa_line{W}(select_tp,:) = collab_loc(1,:); % ��W��Эͬ���µ�λ�ù��ƽ��
            %     [collab_loc,collab_err] = Cal_collab_N(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true,N);  %N�������嶨λ��Ӧ�����
            %     [collab_loc,collab_err] = Cal_collab_logpath(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true);  %��Ӧ�����
            co_err_line(W,select_tp) = collab_err(1); % ÿһ���Ե�Ķ�λ���
        end
        mean_coerr_line(repeat_time,W) = mean(co_err_line(W,:)); %ѡ���W�����Ƶ���Ϊ�ɶ����ʱ��,��λ���
    end
    
    
    %% λ��Эͬ���
    fuse_loc_line = zeros(size(all_loc_test,1),2); % ��λ�㼯�����ʼ��
    for i = 1:size(all_loc_test,1)
        for W = 1:KK
            fuse_loc_line(i,:) = fuse_loc_line(i,:) + aa_line{W}(i,:).*top_similar(W,i);
        end
        fuse_loc_line(i,:) = fuse_loc_line(i,:)/sum(top_similar(1:KK,i)); % ���յĹ��ƽ��
    end
    for i = 1:size(fuse_loc_line,1)
        err_fuse_line(i) = norm(fuse_loc_line(i,:)-all_loc_test(i,:));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% ����ΪԲ�α߽��Эͬ��λ������ж�λ %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% ������֤һ�¶��ڲ��Ե�select_tp�Ķ�λ�����Ƿ�������
    %% ����ÿһ�����Ե����AP-rank������������ԭʼ�ź�ǿ�Ƚ������򣬶��ǿ����ϰ���
    tic;
    [test_ap_rank] = Cal_test_rank(mean_rssi_test);
    
    
    % cc = cc_no;
    %% ȫ�����ж�λ
    for W = 1:KK %ѡ����ٸ�����Эͬ��λ��ǰ��
        for select_tp = 1:size(all_loc_test,1) % ��ǰ��λselect_tp��
            disp(select_tp);
            tp_index = [];
            %     [tp_index,top_similar] = Cal_ap_rank_sim(test_ap_rank,select_tp,mean_rssi_test,N,W);
            [tp_index,top_similar(W,select_tp)] = Cal_ap_rank_sim2(test_ap_rank,select_tp,mean_rssi_test,W);
            tp_index = [select_tp;tp_index]; % �ѳɶԲ��Ե������д����
            collab_rssi = [];
            collab_rssi = mean_rssi_test(tp_index(1),:);
            for i = 2:length(tp_index)
                collab_rssi = [collab_rssi;mean_rssi_test(tp_index(i),:)];% ��Ӧ�ĳɶԲ��Ե��½����ź�ǿ��
            end
            collab_cc = [];
            for i = 1:length(tp_index)
                collab_cc{i} = cc{tp_index(i)};% ����ȷ������cc����cc_no��ÿһ����λ��������λ�ڵ������
            end
            %% �������ҵ���������tp_index��˴β��Ե�����1����peer-to-peer��Эͬ��λ
            tp_collab_true = all_loc_test(tp_index(1),:);
            for i = 2:length(tp_index)
                tp_collab_true = [tp_collab_true;all_loc_test(tp_index(i),:)]; % ��ʵλ��
            end
            
            
            [collab_loc,collab_err] = Cal_collab(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true);  %��Ӧ�����
            aa{W}(select_tp,:) = collab_loc(1,:); % ��W��Эͬ���µ�λ�ù��ƽ��
            %     [collab_loc,collab_err] = Cal_collab_N(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true,N);  %��Ӧ�����
            %     [collab_loc,collab_err] = Cal_collab_logpath(tp_index,collab_cc,vir_rp,collab_rssi,ap,tp_collab_true);  %��Ӧ�����
            co_err(W,select_tp) = collab_err(1); % ÿһ���Ե�Ķ�λ���
            
        end
        mean_coerr(repeat_time,W) = mean(co_err(W,:)); %ѡ���W�����Ƶ���Ϊ�ɶ����ʱ��,Բ�α߽綨λ���
    end
    
    
    %% λ��Эͬ���
    fuse_loc = zeros(size(all_loc_test,1),2); % ��λ�㼯�����ʼ��
    for i = 1:size(all_loc_test,1)
        for W = 1:KK
            fuse_loc(i,:) = fuse_loc(i,:) + aa{W}(i,:).*top_similar(W,i);
        end
        fuse_loc(i,:) = fuse_loc(i,:)/sum(top_similar(1:KK,i)); % ���յĹ��ƽ��
    end
    for i = 1:size(fuse_loc,1)
        err_fuse(i) = norm(fuse_loc(i,:)-all_loc_test(i,:));
    end
    t_cb2 = toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% OCLoc-N��λ���� %%%%%%%%%%%%%%%%%%%%%%%%%
%     tic
%     for select_tp = 1:size(all_loc_test,1)
%         disp(select_tp);
%         [oc_loc_N(select_tp,:),oc_err_N(select_tp)] = OCLoc_N(mean_rssi_no_pen,mean_rssi_test,all_loc,all_loc_test,KK,select_tp,test_ap_rank);
%     end
%     t_ocn = toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% OCLoc-N��λ���� %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% OCLoc-W��λ���� %%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    for select_tp = 1:size(all_loc_test,1)
        disp(select_tp);
        [oc_loc(select_tp,:),oc_err(select_tp)] = OCLoc(mean_rssi_no_pen,mean_rssi_test,all_loc,all_loc_test,KK,select_tp,test_ap_rank);
    end
    t_oc = toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% OCLoc-W��λ���� %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% ����CDFͼ
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
    %     %     lsr{i+4} = ['��',num2str(i)]; % ͼע
    %     % end
    %     KK = 0;
    %     plot(sort(err_fuse_line),[1/125:1/125:1],':k');
    %     lsr{KK+5} = ['ֱ��Эͬ'];
    %     hold on
    %     plot(sort(err_fuse),[1/125:1/125:1],'r');
    %     lsr{KK+6} = ['Բ��Эͬ'];
    %     % legend('CBWF','����','����','����','����','���');
    %     legend(lsr);
    %     xlabel('��λ���ף�');
    %     ylabel('CDF');
    %
    %     %% ������
    %     fprintf('NN�����Ϊ��%.2f\n',mean(est_err_NN));
    %     fprintf('KNN�����Ϊ��%.2f\n',mean(est_err_KNN));
    %     fprintf('ֱ�߽߱�����Ϊ��%.2f��\n',mean(est_err_line));
    %     fprintf('Բ�α߽�����Ϊ��%.2f��\n',mean(est_err));
    %     for i =1:KK
    %         fprintf('��%d��Эͬ���£���Эͬ��λ�����Ϊ��%.2f��\n',i,mean_coerr(i));
    %     end
    %     fprintf('���ֱ��Эͬ��λ�����Ϊ��%.2f��\n',mean(err_fuse_line));
    %     fprintf('���Բ��Эͬ��λ�����Ϊ��%.2f��\n',mean(err_fuse));
    aa_NN(repeat_time) = mean(est_err_NN);
    aa_KNN(repeat_time) = mean(est_err_KNN);
    aa_li(repeat_time) = mean(est_err_line);
    aa_li_coll(repeat_time) = mean(err_fuse_line);
    aa_cb(repeat_time) = mean(est_err);
    aa_cb_coll(repeat_time) = mean(err_fuse);
    % ÿһ��ѭ���µõ���ƽ����λ���
end
%% ����3D������
% figure;
% X = [3.87,3.51,3.75,3.12;3.60,3.11,3.36,2.85;3.53,3.09,3.03,2.62;3.65,3.34,3.02,2.71;3.86,3.78,3.51,3.47]; %(ϵ������)
% E = [0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];%����׼�
% legends = {'Line','Line\_col','CBWF','CBWF\_col'};%��ͼ����
% groupnames = cell(9,1);
% groupnames{1} = '1';groupnames{2} = '2';groupnames{3} = '3';groupnames{4} = '4';%���������ƣ�
% groupnames{5} = '5';
% groupnames{6} = '5';groupnames{7} = '6';groupnames{8} = '7';%���������ƣ�
% Title = '';
% Xlabel = '\tau';
% Ylabel = 'Localization accuracy (m)';
% barweb(X,E,1,groupnames,Title,Xlabel,Ylabel,jet,'none',legends,2,'plot');
% h= legend('Line','Line\_col','CBWF','CBWF\_col');
% box on;
% grid on;
% set(h,'Orientation','horizon');


%% ����CDFͼ
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
% X = [4.97,4.29,4.28;5.02,4.32,4.41;5.08,4.39,4.52;5.09,4.45,4.78;5.44,4.96,5.28;6.50,5.65,6.74;7.50,6.84,7.94;8.32,7.89,9.18]; %(ϵ������)
% E = [0.09,0.07,0.13;0.12,0.12,0.21;0.24,0.18,0.18;0.22,0.23,0.61;0.58,0.38,0.50;0.86,0.49,1.98;1.29,1.25,1.88;1.65,2.03,2.71];%����׼�
% groupnames = cell(8,1);
% groupnames{1} = '0.25';groupnames{2} = '0.5';groupnames{3} = '0.75';groupnames{4} = '1';%���������ƣ�
% groupnames{5} = '2';groupnames{6} = '4';groupnames{7} = '6';groupnames{8} = '8';%���������ƣ�
% Title = '';
% Xlabel = 'SD of Gaussian Noise';
% Ylabel = 'Localization accuracy (m)';
% barweb(X,E,1,groupnames,Title,Xlabel,Ylabel,jet,'none',legends,2,'plot');
% legend('Line','CBWF w/o \omega','CBWF w/ \omega');
% box on;
% grid on;