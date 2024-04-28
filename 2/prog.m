format shortG
for J=[20,40,80,160,320]
    d_x=2*pi/J;    d_t=0.4*d_x*d_x;
    x_array=transpose((-pi+d_x):d_x:pi);
    N=ceil(1/d_t);  T=N*d_t;
    ERR=0.0000001/sqrt(d_x);
    %Advancing Matrixes    niu*a=0.4
    B0=diag(ones(1,J)*0.2)+diag(ones(1,J-1)*0.4,1)+diag(ones(1,J-1)*0.4,-1);
    B0(1,J)=0.4; B0(J,1)=0.4;
    B1=diag(ones(1,J)*1.8)-diag(ones(1,J-1)*0.4,1)-diag(ones(1,J-1)*0.4,-1);
    B1(1,J)=-0.4; B1(J,1)=-0.4;
    BDF=diag(ones(1,J-1)*0.8,1)+diag(ones(1,J-1)*0.8,-1);
    BDF(1,J)=0.8; BDF(J,1)=0.8;
    for m=1:3
        if m==1         %Initial Value 1 + Real Solution
            u0=zeros(J,1);
            for j=1:J
                if(j>=J/4 && j<=3*J/4)
                    u0(j)=1;
                end            
            end
            u1=0.5*ones(J,1);
            u_tmp=ones(J,1);
            i=0;
            while norm(u_tmp)>ERR
                %让级数的相邻误差在10^-7量级，比pde数值解的误差小至少2个量级
                i=i+1; nn=2*i-1;
                u_tmp=2*(2*mod(i,2)-1)*cos(nn*x_array)*exp(-nn*nn*T)/nn/pi;
                u1=u1+u_tmp;
            end
        elseif m==2     %Initial Value 2 + Real Solution
            u0=pi*ones(J,1)-abs(x_array);
            u1=pi*ones(J,1)/2;
            u_tmp=ones(J,1);
            nn=-1;
            while norm(u_tmp)>ERR
                nn=nn+2;
                u_tmp=4*cos(nn*x_array)*exp(-nn*nn*T)/pi/nn/nn;
                u1=u1+u_tmp;
            end
        else            %Initial Value 3 + Real Solution
            u0=cos(x_array);
            u1=exp(-T)*cos(x_array);
        end
        for k=1:2
            u=u0;
            if k==1     %Start with fully explicit scheme
                uu=B0*u;
            else        %Start with fully implicit scheme
                uu=B1\u;
            end
            for n=2:N   %DF scheme (3-Level scheme)
                u=(0.2*u+BDF*uu)/1.8;
                utmp=u; u=uu; uu=utmp;
            end
            [J,m,k,norm(uu-u1)*sqrt(d_x)]
        end
    end
end
