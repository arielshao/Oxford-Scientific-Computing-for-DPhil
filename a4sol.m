function a4sol=a4sol(x1best,y1best,ttbest)
% x1=0.4884; y1=-.0974; tt=3.2504
x1=x1best;y1=y1best; tt=ttbest;
opts=optimset('tolfun',1e-14);
X=fminunc(@dif, [x1;y1;tt],opts);
-dif(X)

function dif=dif(X)
x1=X(1); y1=X(2);tt=X(3);
r1=sqrt(x1^2+y1^2); t1=atan2(y1,x1);
f1=x1.*exp(-r1.^2).*sin(5*(t1+r1));
x2=x1+cos(tt);y2=y1+sin(tt);
r2=sqrt(x2^2+y2^2); t2=atan2(y2,x2);
f2=x2.*exp(-r2.^2).*sin(5*(t2+r2));
dif=-abs(f2-f1);
end
end
    



    