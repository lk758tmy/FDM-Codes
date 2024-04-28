format shortG

J=320;
d_x=2*pi/J;    d_t=0.4*d_x*d_x;
x_array=transpose((-pi+d_x):d_x:pi);
N=ceil(1/d_t); %N=N*10;
ERR=0.0000001/sqrt(d_x);

%Advancing Matrixes    niu*a=0.4
B0=diag(ones(1,J)*0.2)+diag(ones(1,J-1)*0.4,1)+diag(ones(1,J-1)*0.4,-1);
B0(1,J)=0.4; B0(J,1)=0.4;
B1=diag(ones(1,J)*1.8)-diag(ones(1,J-1)*0.4,1)-diag(ones(1,J-1)*0.4,-1);
B1(1,J)=-0.4; B1(J,1)=-0.4;
BDF=diag(ones(1,J-1)*0.8,1)+diag(ones(1,J-1)*0.8,-1);
BDF(1,J)=0.8; BDF(J,1)=0.8;

for m=1:3
    if m==1         %Initial Value 1
        u0=zeros(J,1);
        for j=1:J
            if(j>=J/4 && j<=3*J/4)
                u0(j)=1;
            end            
        end
    elseif m==2     %Initial Value 2
        u0=pi*ones(J,1)-abs(x_array);
    else            %Initial Value 3
        u0=cos(x_array);
    end

    v1=zeros(J+1,1);
    v2=zeros(J+1,1);

    for k=1:2
        u=u0;
        for n=1:N   %DF scheme (3-Level scheme)
            if(n==1)    %Start with fully explicit/implicit
                if k==1
                    uu=B0*u;
                else
                    uu=B1\u;
                end
            else
                u=(0.2*u+BDF*uu)/1.8;
                utmp=u; u=uu; uu=utmp;
            end

            T=n*d_t;
            if m==1     %calculate real solution
                u1=0.5*ones(J,1);
                u_tmp=ones(J,1);
                i=0;
                while norm(u_tmp)>ERR
                    i=i+1; nn=2*i-1;
                    u_tmp=2*(2*mod(i,2)-1)*cos(nn*x_array)*exp(-nn*nn*T)/nn/pi;
                    u1=u1+u_tmp;
                end
            elseif m==2
                u1=pi*ones(J,1)/2;
                u_tmp=ones(J,1);
                nn=-1;
                while norm(u_tmp)>ERR
                    nn=nn+2;
                    u_tmp=4*cos(nn*x_array)*exp(-nn*nn*T)/pi/nn/nn;
                    u1=u1+u_tmp;
                end
            else
                u1=exp(-T)*cos(x_array);
            end

            if k==1     %calculate and record L2 error
                v1(n+1)=norm(uu-u1)*sqrt(d_x);
            else
                v2(n+1)=norm(uu-u1)*sqrt(d_x);
            end
        end
    end

    figure(m)
    semilogy(v1) %plot(v1)
    hold on
    semilogy(v2) %plot(v2)
    title('L2 Error (DF scheme)')
    xlabel('n')
    ylabel('Error')
    legend('Start with fully explicit','Start with fully implicit')
end
