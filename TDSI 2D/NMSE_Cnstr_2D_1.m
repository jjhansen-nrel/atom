function [temp_T,temp_vx,temp_vy]=NMSE_Cnstr_2D_1(pinv_Rdd,R_md,J,C,B,sigmaT2,sigmavx2,sigmavy2)
% NMSE_Cnstr is a service progam for t_stochastic_inverse3;
% it calculates NMSE of the constrained 3D solution.

temp_T=zeros(J,1);
temp_vx=temp_T;
temp_vy=temp_T;
CT=C';
if sigmaT2~=0
    if ~isempty(B) && any(B) 
        for i=1:J
            tmp=CT(i,:)/(C*CT);
            tmp1=tmp*C*R_md(1:J,:);
            tmp2=tmp*B;
            temp_T(i)=1-R_md(i,:)*pinv_Rdd*R_md(i,:)'/sigmaT2+...
                tmp1*pinv_Rdd*tmp1'/sigmaT2+tmp2*tmp2'/sigmaT2;
        end
    else
        for i=1:J
            tmp=CT(i,:)/(C*CT);
            tmp1=tmp*C*R_md(1:J,:);
            temp_T(i)=1-R_md(i,:)*pinv_Rdd*R_md(i,:)'/sigmaT2+...
                tmp1*pinv_Rdd*tmp1'/sigmaT2;
        end
    end
end
if sigmavx2~=0
    if ~isempty(B) && any(B) 
        for i=1:J
            tmp=CT(i,:)/(C*CT);
            tmp1=tmp*C*R_md(J+1:3*J,:);
            tmp2=tmp*B;
            temp_vx(i)=1-R_md(i+J,:)*pinv_Rdd*R_md(i+J,:)'/sigmavx2+...
                tmp1*pinv_Rdd*tmp1'/sigmavx2+tmp2*tmp2'/sigmavx2;
            
        end
    else
        for i=1:J
            tmp=CT(i,:)/(C*CT);
            tmp1=tmp*C*R_md(J+1:3*J,:);
            temp_vx(i)=1-R_md(i+J,:)*pinv_Rdd*R_md(i+J,:)'/sigmavx2+...
                tmp1*pinv_Rdd*tmp1'/sigmavx2;
        end
    end
end
if sigmavy2~=0
     if ~isempty(B) && any(B) 
        for i=1:J
            tmp=CT(i+J,:)/(C*CT);
            tmp1=tmp*C*R_md(J+1:3*J,:);
            tmp2=tmp*B;
            temp_vy(i)=1-R_md(i+2*J,:)*pinv_Rdd*R_md(i+2*J,:)'/sigmavy2+...
                tmp1*pinv_Rdd*tmp1'/sigmavy2+tmp2*tmp2'/sigmavy2;
            
        end
    else
        for i=1:J
            tmp=CT(i+J,:)/(C*CT);
            tmp1=tmp*C*R_md(J+1:3*J,:);
            temp_vy(i)=1-R_md(i+2*J,:)*pinv_Rdd*R_md(i+2*J,:)'/sigmavy2+...
                tmp1*pinv_Rdd*tmp1'/sigmavy2;
        end
    end
end

