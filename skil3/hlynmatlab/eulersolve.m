function k=eulersolve(y0,T,n)
h=T/n;
y=y0;t=0;
yvec=y0;tvec=0;
for i=1:n
    y(:,i+1)=y(:,i)+h*DF(t,y(:,i)) ;
    t=t+h;
    tvec(i+1)=t;
end
k=[tvec ;y];
end