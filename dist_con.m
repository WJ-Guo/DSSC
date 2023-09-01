function [times, is_nei,nei_all] = dist_con(nei, p1, p2, nei_num, max_time)
times = 0;
is_nei = 0;
nei_all = [nei(p1,1:nei_num),p1]; % ��ʼ���ھӼ�
nei_new = nei_all; % �����ж������ھ�
iter = 0;
while is_nei == 0
    if ismember(p2, nei_all)
        is_nei = 1;
        break
    elseif isempty(nei_new)
        break
    end
    times = times+1;
    
    nei_temp = nei_new;
    nei_new = [];
    dim = size(nei_temp, 2);
    % ��ȡÿһ�����ڵĽ���
    for i = 1:dim
       p3 = nei_temp(i);
       nei_p3 = nei(p3,:);
       for j = 1:length(nei_p3)
           if ~ismember(nei_p3(j), nei_all)
              nei_new = [nei_new, nei_p3(j)];
              nei_all = [nei_all, nei_p3(j)];
           end
       end
    end
    
    if times >  max_time
        break
    end
    iter = iter +1;
end


