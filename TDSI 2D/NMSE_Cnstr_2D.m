function [temp_T,temp_vx,temp_vy]=NMSE_Cnstr_2D(G_1,R_md,J,D1,D2,D3,F1,F2,F3,sigmaT2,sigmavx2,sigmavy2)
% NMSE_Cnstr is a service progam for t_stochastic_inverse3;
% it calculates NMSE of the constrained 3D solution.

temp_T=zeros(J,1);
temp_vx=temp_T;
temp_vy=temp_T;
for j1=1:3
    for j2=1:3
        GR=G_1((j1-1)*J+1:j1*J,:)*R_md((j2-1)*J+1:j2*J,:)';
        if j1==j2 && j1==1
            if sigmaT2~=0
                diag_T=diag(GR);
                temp_T=temp_T+diag(D1*GR*D1');
            end
        end
        if j1==j2 && j1==2
            diag_vx=diag(GR);
        end
        if j1==j2 && j1==3
            diag_vy=diag(GR);
        end
        if sigmavx2~=0 && j1>1 && j2>1
            temp_vx=temp_vx+diag(D2(:,(j1-2)*J+1:(j1-1)*J)*GR*D2(:,(j2-2)*J+1:(j2-1)*J)');
        end
        if sigmavy2~=0 && j1>1 && j2>1
            temp_vy=temp_vy+diag(D3(:,(j1-2)*J+1:(j1-1)*J)*GR*D3(:,(j2-2)*J+1:(j2-1)*J)');
        end
    end
end
clear D1 D2 D3 R_md G_1 GR;
if sigmaT2~=0
    if isempty(F1)
        temp_T=ones(J,1)+(temp_T-diag_T)/sigmaT2;
    else
        temp_T=ones(J,1)+(temp_T-diag_T+F1.^2)/sigmaT2;
    end
end
if sigmavx2~=0
    if isempty(F2)
        temp_vx=ones(J,1)+(temp_vx-diag_vx)/sigmavx2;
    else
        temp_vx=ones(J,1)+(temp_vx-diag_vx+F2.^2)/sigmavx2;
    end
end
if sigmavy2~=0
    if isempty(F3)
        temp_vy=ones(J,1)+(temp_vy-diag_vy)/sigmavy2;
    else
        temp_vy=ones(J,1)+(temp_vy-diag_vy+F3.^2)/sigmavy2;
    end
end

