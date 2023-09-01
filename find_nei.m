function [nei_all] = find_nei(nei, p1, nei_num, max_time)
times = 0;
nei_all{1} = [nei(p1,1:nei_num),p1]; % 初始化邻居集
nei_new = nei_all{1}; % 用于判断新增邻居

while true    

    times = times+1;
    if times >= max_time
        break
    end
    
    nei_temp = nei_new;
    nei_new = [];
    dim = size(nei_temp, 2);
    nei_all{times+1} = nei_all{times};
    % 读取每一个近邻的近邻
    for i = 1:dim
       p3 = nei_temp(i);
       nei_p3 = nei(p3,1:nei_num);
       for j = 1:length(nei_p3)
           if ~ismember(nei_p3(j), nei_all{times+1})
              nei_new = [nei_new, nei_p3(j)];
              nei_all{times+1}= [nei_all{times+1}, nei_p3(j)];
           end
       end
    end
     
end


