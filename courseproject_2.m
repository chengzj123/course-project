tic
a=-1;
c=1;
Tol=10^(-4);
NN=5;
r=10^(-NN);
C=24;  
f=@(x) 0;
p=@(x) 2*x/r;
q=@(x) 0;
GL=@(x) x+1;
GR=@(x) x-1;
s=2;
PhiL=@(x) x/r;
PhiR=@(x) x/r;
ui=@(x) x;
F=@(x) -2*x/r;
[Ur,eI]=linearsolver(a,c,Tol,C,GL,GR,s,PhiL,PhiR,ui,F,r);
plot(Ur(1,:),Ur(2,:),'k');
toc
disp( ['运行时间: ',num2str(toc) ] );




