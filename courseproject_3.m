tic
Fun=struct('U',@(x) 0);
dFun=struct('V',@(x) 0);
ddFun=struct('R',@(x) 0);
Fun(1)=struct('U',@(x) (x+2*exp(2)-2*exp(1))/(2*exp(2)-exp(1)));
dFun(1)=struct('V',@(x) 1/(2*exp(2)-exp(1)));
ddFun(1)=struct('R',@(x) 0);
f=@(x,y,z) z/(x*(y+1)^2)-y/(x^2*(y+1));
f2=@(x,y,z) -y/(x^2*(y+1)^2)-2*z/(x*(y+1)^3);
f3=@(x,y,z) 1/(x*(y+1)^2);
X=exp(1):(2*exp(2)-exp(1))/100:2*exp(2);
Y=X;
tol=10^(-4);
er=1;
N=1;
a=exp(1);
c=2*exp(2);
C=16;
while er>tol
      e=0;
      f1=0;
      p=@(x) -f3(x,Fun(N).U(x),dFun(N).V(x));
      q=@(x) -f2(x,Fun(N).U(x),dFun(N).V(x));
      F=@(x) -ddFun(N).R(x)+f(x,Fun(N).U(x),dFun(N).V(x));
      GL=@(x) x-exp(1);
      GR=@(x) x-2*exp(2);
      dGL=@(x) 1;
      dGR=@(x) 1;
      s=2*exp(2)-exp(1);
      PhiL=@(x) (p(x)+q(x)*GR(x))/s;
      PhiR=@(x) (p(x)+q(x)*GL(x))/s;
      [Ur,JL,JR,Leftendpoint]=linearsolver1(a,c,tol,C,GL,GR,s,PhiL,PhiR,F);
      d2v=@(x) f3(x,Fun(N).U(x),dFun(N).V(x))*dv(x,Ur,JL,JR,Leftendpoint...
          ,GL,GR,dGL,dGR,s,c)+ f2(x,Fun(N).U(x),dFun(N).V(x))*v(x,Ur,JL,...
          JR,Leftendpoint,GL,GR,s,c)+F(x);
      Fun(N+1)=struct('U',@(x) 0);
      Fun(N+1).U=@(x) Fun(N).U(x)+v(x,Ur,JL,JR,Leftendpoint,GL,GR,s,c);
      dFun(N+1)=struct('V',@(x) 0);
      dFun(N+1).V=@(x) dFun(N).V(x)+dv(x,Ur,JL,JR,Leftendpoint,GL,GR,...
          dGL,dGR,s,c);
      ddFun(N+1)=struct('R',@(x) 0);
      ddFun(N+1).R=@(x) ddFun(N).R(x)+d2v(x);      
      for i=1:101
          e=max(e,abs(Fun(N+1).U(X(i))-Fun(N).U(X(i))));
          f1=max(f1,abs(Fun(N+1).U(X(i))+Fun(N).U(X(i))));
      end
      er=e/f1;
      N=N+1;
      if N==5
         break
      end
end
for j=1:4
    for i=1:101
        Y(j,i)=Fun(j).U(X(i));
    end
end
for i=1:101
    Y(3,i)=lambertw(X(i));
end
ER=0;
for i=1:101
    ER=max(ER,abs(Y(4,i)-lambertw(X(i))));
end
plot (X,Y(1,:),'k',X,Y(2,:),'r',X,Y(3,:),'g')
toc
disp( ['运行时间: ',num2str(toc) ] );

function deltau=v(x,Ur,JL,JR,Leftendpoint,GL,GR,s,c)
         n=length(Leftendpoint)+1;
         YY=zeros(1,n);
         for i=1:n-1
             YY(i)=Leftendpoint(i);
         end
         YY(n)=c;
         for i=1:n-1
             if x>=YY(i) && YY(i+1)>=x
                K=Ch([YY(i),YY(i+1)]);
                Z1=zeros(10,1);
                Z2=zeros(10,1);
                for j=1:10
                    Z1(j)=GL(K(j));
                    Z2(j)=GR(K(j));
                end
                deltau=GR(x)/s*(JL(i)+LIntegral(YY(i),x,YY(i+1),Z1.*...
                    (Ur(2,10*(i-1)+1:10*i))'))+GL(x)/s*(JR(i+1)+ ...
                       RIntegral(YY(i),x,YY(i+1),Z2.*(Ur(2,10*(i-1)+...
                       1:10*i))'));
                break
             end
         end
end

function deltadu=dv(x,Ur,JL,JR,Leftendpoint,GL,GR,dGL,dGR,s,c)
         n=length(Leftendpoint)+1;
         YY=zeros(1,n);
         for i=1:n-1
             YY(i)=Leftendpoint(i);
         end
         YY(n)=c;
         for i=1:n-1
             if x>=YY(i) && YY(i+1)>=x
                K=Ch([YY(i),YY(i+1)]);
                Z1=zeros(10,1);
                Z2=zeros(10,1);
                for j=1:10
                    Z1(j)=GL(K(j));
                    Z2(j)=GR(K(j));
                end
                deltadu=dGR(x)/s*(JL(i)+LIntegral(YY(i),x,YY(i+1),Z1.*...
                    (Ur(2,10*(i-1)+1:10*i))'))+dGL(x)/s*(JR(i+1)+ ...
                       RIntegral(YY(i),x,YY(i+1),Z2.*(Ur(2,10*(i-1)+...
                       1:10*i))'));
                break
             end
         end
end


         






