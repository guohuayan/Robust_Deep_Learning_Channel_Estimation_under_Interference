function local = rand_ue_location(num,center,R)
    xy=zeros(num,3);    
    for i=1:num
        [xy(i,1),xy(i,2)]=RandCircle(R);     %注意要分开赋值！！
    end
    local=xy+ones(num,1)*center;
end

%% 反函数法
function [x,y]=RandCircle(radius)
    rand_theta=rand()*pi*2;
    rand_radius=sqrt(rand())*radius;
    x=cos(rand_theta)*rand_radius;
    y=sin(rand_theta)*rand_radius;
end

%% 舍选法
% function [x,y]=RandCircle(radius)
%     while 1
%         x=(rand()-0.5)*2*radius;
%         y=(rand()-0.5)*2*radius;
%         if (x^2+y^2)<radius^2
%             break;
%         end
%     end
% end