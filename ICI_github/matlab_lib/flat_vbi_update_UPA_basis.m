function [At,beta1,beta2,flag] = flat_vbi_update_UPA_basis(y1,y2,he1,he2,C_old,a1,a2,beta1,beta2,xo1,xo2,At)
    M1=length(beta1);
    M2=length(beta2);
    Ag=he1*he1'+C_old*2+he2*he2';
%     Ag=C_old*2;
%     Ag=he1*he1'+he2*he2';
    Bg=he1*y1'+he2*y2';
%     Bg=0;
%     Ag=-Ag;
    Ai=At;
    f0=real(trace(Ag*(Ai'*Ai)-2*Bg*Ai));
%     f0=real(trace(alpha0*Ag*(Ai'*Ai)-2*alpha0*Bg*Ai+Ag*Ai'*Gh_old*Ai));
    x1=xo1+beta1;
    x2=xo2+beta2;
    beta1_m=beta1;
    beta2_m=beta2;
    flag=0;
    min_beta1=-1-xo1+eps;
    max_beta1=1-xo1-eps;
    min_beta2=-1-xo2+eps;
    max_beta2=1-xo2-eps;
    step=1e-3/max(M1,M2);
    while(1)
        grad_g=2*Ai*Ag-2*Bg';
%         grad_g=2*alpha0*Ai*Ag-2*alpha0*Bg'+2*Gh_old*Ai*Ag;
        A1=exp(a1*x1.')./sqrt(M1);
        A2=exp(a2*x2.')./sqrt(M2);
        A1g=A1.*(a1*ones(1,M1));
        A2g=A2.*(a2*ones(1,M2));
        g1=zeros(M1,1);
        for gx0=1:M1
            ag=zeros(M1,M1);
            ag(:,gx0)=A1g(:,gx0);
            g1(gx0)=real(trace(grad_g'*kron(ag,A2)));
        end
        g2=zeros(M2,1);
        for gx0=1:M2
            ag=zeros(M2,M2);
            ag(:,gx0)=A2g(:,gx0);
            g2(gx0)=real(trace(grad_g'*kron(A1,ag)));
        end
        gc=[g1;g2];
        gc=gc./norm(gc,2);
        g1=gc(1:M1);
        g2=gc((M1+1):end);
%         step=1e-3/max(M1,M2);
        beta1_m=beta1-(g1).*step;
        beta1_m=min(beta1_m,1/M1);beta1_m=max(beta1_m,-1/M1);
%         beta1=min(beta1,max_beta1);beta1=max(beta1,min_beta1);
        beta2_m=beta2-(g2).*step;
        beta2_m=min(beta2_m,1/M2);beta2_m=max(beta2_m,-1/M2);
%         beta2=min(beta2,max_beta2);beta2=max(beta2,min_beta2);
        x1=xo1+beta1_m;
        x2=xo2+beta2_m;
        A1=exp(a1*x1.')./sqrt(M1);
        A2=exp(a2*x2.')./sqrt(M2);
        Ai=kron(A1,A2);
        f1=real(trace(Ag*(Ai'*Ai)-2*Bg*Ai));
%         (f1-f0)/abs(f0)
%         pause
%         f1=real(trace(alpha0*Ag*(Ai'*Ai)-2*alpha0*Bg*Ai+Ag*Ai'*Gh_old*Ai));
        if f1>f0 
            if flag==0
                step=step/2;
%                 step
            else
                x1=xo1+beta1;
                x2=xo2+beta2;
                A1=exp(a1*x1.')./sqrt(M1);
                A2=exp(a2*x2.')./sqrt(M2);
                Ai=kron(A1,A2);
                At=Ai;
                break
            end
        elseif flag==1 && abs(f1-f0)/abs(f0)<1e-6
            beta1=beta1_m;
            beta2=beta2_m;
            x1=xo1+beta1;
            x2=xo2+beta2;
            A1=exp(a1*x1.')./sqrt(M1);
            A2=exp(a2*x2.')./sqrt(M2);
            Ai=kron(A1,A2);
            At=Ai;
            break
        else
            beta1=beta1_m;
            beta2=beta2_m;
            f0=f1;
            flag=1;
        end
    end     
end

