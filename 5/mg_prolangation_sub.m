function I=mg_prolangation_sub(height,width,p)
I=sparse(height-width,width);
for i=1:size(p,2)
    pp=p(:,i);
    p1=pp(1); p2=pp(2); p3=pp(3);
    if(p1>width)
        if(p2<=width)
            I(p1-width,p2)=0.5;
        end
        if(p3<=width)
            I(p1-width,p3)=0.5;
        end
    end
    if(p2>width)
        if(p1<=width)
            I(p2-width,p1)=0.5;
        end
        if(p3<=width)
            I(p2-width,p3)=0.5;
        end
    end
    if(p3>width)
        if(p1<=width)
            I(p3-width,p1)=0.5;
        end
        if(p2<=width)
            I(p3-width,p2)=0.5;
        end
    end
end

