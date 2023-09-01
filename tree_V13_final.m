clear all; 
clc; 
close all;
load A2.txt;
data = A2(:,1:(end-1));
tic
data = data - repmat(mean(data),size(data,1),1);
data = data/max(max(abs(data)));
			
true_label = A2(:,end);


NClusters = max(true_label);

figure
hold on
cmap=colormap;
for i=1:NClusters
       ic=int8((i*64.)/(NClusters*1.));

       if NClusters<=12
               switch i
                case 1
                      plot(data(true_label==i,1),data(true_label==i,2),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 2
                      plot(data(true_label==i,1),data(true_label==i,2),'*','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 3
                      plot(data(true_label==i,1),data(true_label==i,2),'x','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 4
                      plot(data(true_label==i,1),data(true_label==i,2),'+','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 5
                      plot(data(true_label==i,1),data(true_label==i,2),'s','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 6
                      plot(data(true_label==i,1),data(true_label==i,2),'d','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 7
                      plot(data(true_label==i,1),data(true_label==i,2),'v','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 8
                      plot(data(true_label==i,1),data(true_label==i,2),'^','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 9
                      plot(data(true_label==i,1),data(true_label==i,2),'<','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 10
                      plot(data(true_label==i,1),data(true_label==i,2),'>','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 11
                      plot(data(true_label==i,1),data(true_label==i,2),'p','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 12
                      plot(data(true_label==i,1),data(true_label==i,2),'h','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
               end
               
       else
            plot(data(true_label==i,1),data(true_label==i,2),'.','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));	
       end % end of if 					  

end %end of for
title('A2','FontSize',12.0)	
hold off

% build kd-tree
NS = createns(data,'NSMethod','kdtree');

dim = size(data,2);
num = size(data,1);

% search knn for every points
knn = 50;
find_k = 16;
max_time = 2; % �����ж��Ƿ�Ϊ�ڽ�Ⱥ
max_nei = 10;
num_k = knn + 1 ;  %knn+1,for thr first con
nei = zeros(num, num_k); % record nei index
dist = zeros(num, num_k); % record nei distance
min_ther = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����ɸѡ��ѡ���ĵ�
% ��������ÿ������k�����Լ�k���������ƽ�����룬�����ܶ�
for i = 1:num
    [nei(i,:), dist(i,:)] = knnsearch(NS,data(i,:),'k',num_k);
    mean_dist(i,:) = mean(dist(i,:));
end
nei(:,1) = []; % �����õ㱾��
dist(:,1) = [];

% ���ܶ����򣬴��ܶȴ�����±�����ȷ���ֲ�����
[mean_dist_sort, id_all_sort] = sort(mean_dist, 'ascend'); % ���е����ƽ����������

center_id1 = [];
for i = 1: num 
    p1 = id_all_sort(i);
    dist_p1 = mean_dist_sort(i);
    same_point = intersect(nei(p1,1:find_k),center_id1);
    dist_nei = mean_dist(nei(p1,1:find_k));
   
     % ����õ����󲿷��ھӵ���ܶȶ�С������Ϊcan
    if isempty (same_point) && (dist_p1 <= mean(dist_nei))
        center_id1 = [center_id1 ,p1];
    end
    
end   

fprintf('step 1 down: knn \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
center_can = center_id1;
num_can = length(center_can);
fprintf('step 2 down: reverse select\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ѡȡ����
% ��¼ÿһ��can��N����Ѱ����
% ������can��һһѰ���ص�����
% �������ʱ�ޱ�ǣ��򽫸õ���Ϊ���ģ����б�ǣ�����Ϊ������

for i= 1:num_can    
    p1 = center_can(i);
    [nei_temp] = find_nei(nei, p1, max_nei, max_time);
    nei_all{i} = nei_temp{2};% ��¼ÿһ��can��N����Ѱ����
    nei_1(i,:) = [p1,nei(p1,:)];
end

cl_can = zeros(num_can, 1);
center_id = 0;
center_temp = [];
cl = zeros(num,1) -1; % ��¼����center���ڵ���
center_no = [];

for i= 1:num_can
   
    p1 = i; % ���ܶ����Ŀ�ʼѰ��
    dist_p1 = mean_dist(center_can(p1));
    if cl_can(p1) == 0
       center_id = center_id + 1;
       cl_can(p1) = center_id;
       center_temp = [center_temp; center_can(p1)];
    end

    for j = i+1:num_can
        p2 = j;
        dist_p2 = mean_dist(center_can(p2));
             
        nei_same2 = intersect(nei_all{p1},nei_all{p2});
 
        if isempty(nei_same2)
            continue
        end
        
        if cl_can(p2) == cl_can(p1)
            continue
        end
        
        dist_same = mean_dist(nei_same2);

        mean_dist_same = mean(dist_same);
        min_diast_same = min(dist_same);
        max_dist_center = max(dist_p2, dist_p1);  
        mean_dist_center = (dist_p2+dist_p1)/2;
        num1 = find(dist_same<max_dist_center);
        ther = (length(nei_all{p1}) + length(nei_all{p2})) / 5;
        
        length_same = length(nei_same2);
        
        turn = 0;
        if length_same >= ther
            turn = 1;
        elseif length_same > min_ther && mean_dist_same < max_dist_center
            turn = 1;
        elseif length_same > ther/2 && min_diast_same < max_dist_center
            turn = 1;
        end
        
        
        if turn == 1
     
           if cl_can(p2) == 0
              cl_can(p2) = cl_can(p1);
                            
           else
		   
              if  cl_can(p1) > cl_can(p2)
                  point_p1 = find(cl_can == cl_can(p1));
                  cl_can(point_p1) = cl_can(p2);
                  center_p1 = center_can(point_p1);
                  center_temp = setdiff(center_temp, center_p1);
              else
                  point_p2 = find(cl_can == cl_can(p2));
                  cl_can(point_p2) = cl_can(p1);
                  center_p2 = center_can(point_p2);
                  center_temp = setdiff(center_temp, center_p2);
              end
			  
           end
        
        end
        
    end
    
end

center_is = center_temp;
for i = 1:num_can
    cl(center_can(i)) = cl_can(i);
end

fprintf('step 3 down: select center can\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% ����ǩ��������1234
NClusters = max(cl);
id = [];
ture_length = length(unique(cl_can));
for i = 1:ture_length
    index = find(cl == i);
    if isempty(index)
        id = [id,i];
    end
end

j = 1; 
for i = (ture_length + 1) : NClusters
    index = find(cl == i);
    if ~isempty(index)
        cl(index) = id(j);
        j = j + 1;
    end
end

NClusters = max(cl);
fprintf('step 4 down: select real center\n');

%% ���Ǽܵľ���ͼ
label = cl;
figure
hold on
cmap=colormap;
for i=1:NClusters
       ic=int8((i*64.)/(NClusters*1.));

       if NClusters<=12
               switch i
                case 1
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 2
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 3
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 4
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 5
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 6
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 7
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 8
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 9
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 10
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 11
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 12
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                end

       else
            plot(data(label==i,1),data(label==i,2),'.','MarkerSize',6,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));	
       end % end of if 			  

end %end of for
hold on
id = center_is(:,1);
for i = 1:length(id)
%     i
    plot(data(id(i),1),data(id(i),2), '*','MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','k');
    plot(data(id(i),1),data(id(i),2), 'o','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');
%     pause
end
title('A2','FontSize',18.0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����ʣ���

% can��1�������ֱ�ӷ����can��ʣ��㰴�������ԭ�����
for i = 1:num_can
    nei_i = [];
    p1 = center_can(i);
    cl_p1 = cl(p1);
    nei_i = nei(p1,1:max_nei);
    cl(nei_i) = cl_p1;
end
nei_num = 0;

% ��can��һ������ڣ���������������can����������Ĳ��ֵ���Կ�������Ϊ����
halo_can = find(cl == -1);
data_can = data(center_can,:);
for i = 1:length(halo_can)
    p1 = halo_can(i);
    data_p1 = data(p1,:);    
    for j = 1:num_can
        dist_temp(j) = sqrt(sum((data_p1 - data_can(j,:)) .^ 2));
    end
    [dist_p1(i), id] = min(dist_temp); 
    cl(p1) = cl(center_can(id));
end
fprintf('step 5 down: clustering');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��ͼ
label = cl;
figure
hold on
cmap=colormap;
for i=1:NClusters
       ic=int8((i*64.)/(NClusters*1.));

       if NClusters<=12
               switch i
                case 1
                      plot(data(label==i,1),data(label==i,2),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 2
                      plot(data(label==i,1),data(label==i,2),'*','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 3
                      plot(data(label==i,1),data(label==i,2),'x','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 4
                      plot(data(label==i,1),data(label==i,2),'+','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 5
                      plot(data(label==i,1),data(label==i,2),'s','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 6
                      plot(data(label==i,1),data(label==i,2),'d','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 7
                      plot(data(label==i,1),data(label==i,2),'v','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 8
                      plot(data(label==i,1),data(label==i,2),'^','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 9
                      plot(data(label==i,1),data(label==i,2),'<','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 10
                      plot(data(label==i,1),data(label==i,2),'>','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 11
                      plot(data(label==i,1),data(label==i,2),'p','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                case 12
                      plot(data(label==i,1),data(label==i,2),'h','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
                end

       else
            plot(data(label==i,1),data(label==i,2),'.','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));	
       end % end of if 			  

end %end of for

hold on
id = center_is(:,1);
for i = 1:length(id)
%     i
    plot(data(id(i),1),data(id(i),2), '*','MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor','k');
    plot(data(id(i),1),data(id(i),2), 'o','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');
%     pause
end

title('A2','FontSize',18.0)
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FMI,ARI,NMI  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DPCFMI = FMI(true_label(:),cl(:));
fprintf('FMI value of DPC on A2 dataset %i \n', DPCFMI);
[cluster_acc,randindex,DPCARI] = ARI(true_label(:),cl(:));
fprintf('ARI value of DPC on A2 dataset %i \n', DPCARI);
DPCNMI = NMI(true_label(:),cl(:));
fprintf('NMI value of DPC on A2 dataset %i \n', DPCNMI);