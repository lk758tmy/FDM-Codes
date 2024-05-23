J0=500; T1=0.25; T2=1;
for nu=[0.5,0.8] %0.7,0.9]
    nunu=nu*nu; nudiv=(1-nu)/(1+nu);
    d_x=1/J0; d_t=nu*d_x;
    N1=round(T1/d_t); N2=round(T2/d_t); J=round(J0+1+nu*N2);
    u0=zeros(J+1,1);
    for i=1:J
        x=(i-1)*d_x;
        if x<0.1
            u0(i,1)=0;
        elseif x<0.26
            u0(i,1)=(x-0.1)/0.16;
        elseif x<0.42
            u0(i,1)=(0.34-x)/0.08;
        elseif x<0.58
            u0(i,1)=(0.5-x)/0.08;
        elseif x<0.74
            u0(i,1)=-1;
        elseif x<0.9
            u0(i,1)=1;
        else
            u0(i,1)=0;
        end
    end
    for m=1:5
        figure(round(m+(nu-0.5)*50/3))
        plot(u0(1:J0,1))
        hold on
        u=zeros(J+1,1); utmp=u0;
        if m==5
            utmp2=zeros(J+1,1);
            for k=1:J0
                utmp2(k,1)=u0(k+1,1);
            end
        end
        for n=1:N2
            switch m
                case 1 %Upwind
                    u(1,1)=0; %u(1,1)=utmp(1,1);
                    for k=2:min((J0+n),J)
                        u(k,1)=(1-nu)*utmp(k,1)+nu*utmp(k-1);
                    end
                case 2 %Lax-Wendroff
                    u(1,1)=0;
                    %u(1,1)=(1-nunu)*utmp(1,1)+0.5*(nunu-nu)*utmp(2,1);
                    for k=2:min((J0+n),J)
                        u(k,1)=(1-nunu)*utmp(k,1)+0.5*((nunu+nu)* ...
                            utmp(k-1,1)+(nunu-nu)*utmp(k+1,1));
                    end
                case 3 %Lax-Friedrichs
                    u(1,1)=0; %u(1,1)=0.5*(1-nu)*utmp(2,1);
                    for k=2:min((J0+n),J)
                        u(k,1)=0.5*((1+nu)*utmp(k-1,1)+(1-nu)*utmp(k+1,1));
                    end
                case 4 %Thomee (Box)
                    u(1,1)=0;
                    for k=2:min((J0+n),J)
                        u(k,1)=utmp(k-1,1)+nudiv*(utmp(k,1)-u(k-1,1));
                    end
                case 5 %FrogLeap
                    u(1,1)=0; %u(1,1)=utmp2(1,1)-nu*utmp(2,1);
                   for k=2:min((J0+n),J)
                        u(k,1)=utmp2(k,1)+nu*(utmp(k-1,1)-utmp(k+1,1));
                   end
            end
            if n==N1
                plot(u((1+round(nu*N1)):(J0+round(nu*N1)),1))
                hold on
            end
            if m==5
                    utmp2=utmp;
            end
            utmp=u;
        end
        plot(u((1+round(nu*N2)):(J0+round(nu*N2)),1))
        legend('t=0','t=0.25s','t=1s','Location','SouthEast')
        switch m
            case 1
                cmt='UpW ';
            case 2
                cmt='LW ';
            case 3
                cmt='Lax ';
            case 4
                cmt='Box ' ;
            case 5
                cmt='Frog ';
        end
        if nu<0.6
            cmt=[cmt 'nu=0.5'];
        %elseif nu<0.8
        %    cmt=[cmt 'nu=0.7'];
        else
            cmt=[cmt 'nu=0.8']; %cmt=[cmt 'nu=0.9'];
        end
        text(50,-0.5,cmt,'FontSize',13)
        text(50,-0.75,'delta_x=1/500','FontSize',13)
        print(m,'-djpeg',[cmt '.jpg'])
    end
end
