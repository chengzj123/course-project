tic
a=0;
c=2;
Tol=10^(-4);
C=16;  
f=@(x) 2;
p=@(x) -x;
q=@(x) 2;
GL=@(x) x;
GR=@(x) x-2;
s=2;
PhiL=@(x) (x-4)/2;
PhiR=@(x) x/2;
ui=@(x) 2*x;
F=@(x) 2-2*x;
A=zeros(1,101);
for i=1:length(A)
    A(i)=(i-1)*2/(length(A)-1);
end
K=10;
Ur=directsolver(A,GL,GR,s,PhiL,PhiR,ui,F,K);
plot(Ur(1,:),Ur(2,:),'k');
er=0;
for i=1:length(Ur)
    er=max(er,abs(Ur(2,i)-Ur(1,i)^2));
end
toc
disp( ['运行时间: ',num2str(toc) ] );