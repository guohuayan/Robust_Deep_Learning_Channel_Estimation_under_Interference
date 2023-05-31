function [x,y] = random_user_hexagon_drange(r,theta_range,r_range)
    %% r is the radius of cell
    %% theta_range=[theta_min,theta_max] is the range of splitting in [0,360]
    N_point=1e2;
%     if theta_range(1)<0
%         theta_range(1)=theta_range(1)+360;
%     end
%     if theta_range(2)<0
%         theta_range(2)=theta_range(2)+360;
%     end
    while(1)
        xy_temp=rand(2,N_point).*2-1;
        xy_temp=xy_temp.*r;
        [theta,rho] = cart2pol(xy_temp(1,:),xy_temp(2,:));
        idx=find(rho>=r);
        theta(idx)=[];
        rho(idx)=[];
        idx=find(rho>=r_range(2));
        theta(idx)=[];
        rho(idx)=[];
        idx=find(rho<=r_range(1));
        theta(idx)=[];
        rho(idx)=[];
        len=length(theta);
        if len==0
            continue
        end
        %%
        thetad=rad2deg(theta);
        idx=find(thetad<0);
        thetad(idx)=thetad(idx)+360;
        if theta_range(1)<theta_range(2)
            idx1=find(thetad<theta_range(1));
            idx2=find(thetad>theta_range(2));
            idx=union(idx1,idx2);
            theta(idx)=[];
            rho(idx)=[];
            thetad(idx)=[];
        else
            idx1=find(thetad>theta_range(1));
            idx2=find(thetad<theta_range(2));
            idx=union(idx1,idx2);
            theta=theta(idx);
            rho=rho(idx);
            thetad=thetad(idx);
        end
        len=length(theta);
        if len==0
            continue
        end
        %%
        ruler_line_p=0:60:360;        
%         thetad=rad2deg(theta);
        %%
        if len==0
            continue
        else
            diff_thetad=zeros(1,len);
            for i0=1:len
                tmp=thetad(i0);
%                 if tmp>0
                    tmp_out=tmp-ruler_line_p;
%                 else
%                     tmp_out=tmp+ruler_line_p;
%                 end
                diff_thetad(i0)=min(abs(tmp_out));
            end
            r_thr=r./cosd(diff_thetad).*sqrt(3/4);
            idx=find(rho>=r_thr);
            theta(idx)=[];
            rho(idx)=[];
            len=length(theta);
            if len==0
                continue
            else
                break
            end
        end
    end
    [xw,yw] = pol2cart(theta,rho);
    x=xw(1);
    y=yw(1);
end

% function flag=cross_2D(x,y)
%     z=x(1)*y(2)-x(2)*y(1);
%     flag=sign(z);
% end