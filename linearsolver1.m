function [Ur,JL,JR,Leftendpoint]=linearsolver1(a,c,Tol,C,GL,GR,s,PhiL,PhiR,F)
         interval=struct('leftchild',0,'rightchild',0,'parent',0,...
             'content',0,'exist',1);
         interval.content=struct('Interval',[0,0],'alphaL',0,'alphaR',...
             0,'betaL',0,'betaR',0, 'deltaL',0,'deltaR',0,'niuL',0, ...
         'niuR',0,'niu',0);
         interval(1).content.Interval=[a,c]; 
         interval(1).content.niuL=0;
         interval(1).content.niuR=0;  
         interval(1).content.niu=1;
         N=1;
         Er=1;        
         Mat2=zeros(4,3);
         Mat2(1,3)=0.5;
         while Er>Tol
               S=zeros(1,N);
               for i=1:N
                   if interval(i).leftchild==0 && interval(i).exist==1                   
                      S(i)=Evamonitor(interval,i,PhiL,PhiR,GL,GR,F);
                   end
               end
               M=N;
               for i=1:N
                   if interval(i).leftchild==0 && interval(i).exist==1
                      if S(i)>=max(S)/C
                         x=interval(i).content.Interval(1);
                         y=interval(i).content.Interval(2);
                         M=M+1;
                         interval(M)=struct('leftchild',0,'rightchild',...
                             0,'parent',i,'content',0,'exist',1);
                         interval(M).content=struct('Interval',[x,...
                             (x+y)/2],'alphaL',0, ...
                         'alphaR',0,'betaL',0,'betaR',0, 'deltaL',0,...
                         'deltaR',0,'niuL',0, ...
                         'niuR',0,'niu',0);               
                         interval(i).leftchild=M;
                         M=M+1;
                         interval(M)=struct('leftchild',0,'rightchild',...
                             0,'parent',i,'content',0,'exist',1);
                         interval(M).content=struct('Interval',[(x+y)/2,...
                             y],'alphaL',0, ...
                         'alphaR',0,'betaL',0,'betaR',0, 'deltaL',0,...
                         'deltaR',0,'niuL',0,'niuR',0,'niu',0);
                         interval(i).rightchild=M;            
                      elseif interval(i).parent>0
                             n=interval(i).parent;
                             m=interval(n).leftchild;
                             if S(m)+S(m+1)<max(S)/(2^10) && ...
                                 interval(2*m+1-i).leftchild==0
                                interval(m).exist=0;
                                interval(m+1).exist=0;
                                interval(n).leftchild=0;
                                interval(n).rightchild=0;
                             end                 
                       end
                   end
               end     
               for i=1:M
                   if interval(i).leftchild==0 && interval(i).exist==1           
                      X=ABD(interval,i,F,PhiL,PhiR,GL,GR);
                      interval(i).content.alphaL=X(1);
                      interval(i).content.alphaR=X(2);
                      interval(i).content.betaL=X(3);
                      interval(i).content.betaR=X(4);
                      interval(i).content.deltaL=X(5);
                      interval(i).content.deltaR=X(6); 
                   end
               end
               for j=1:M
                   i=M+1-j;
                   n=interval(i).parent;
                   if n>0 && interval(i).exist==1
                      m=interval(n).leftchild;
                      U(1)=interval(m).content.alphaL;
                      U(2)=interval(m).content.alphaR;
                      U(3)=interval(m).content.betaL;
                      U(4)=interval(m).content.betaR;
                      U(5)=interval(m).content.deltaL;
                      U(6)=interval(m).content.deltaR;
                      V(1)=interval(m+1).content.alphaL;
                      V(2)=interval(m+1).content.alphaR;
                      V(3)=interval(m+1).content.betaL;
                      V(4)=interval(m+1).content.betaR;
                      V(5)=interval(m+1).content.deltaL;
                      V(6)=interval(m+1).content.deltaR;
                      Delta=1-U(3)*V(2);
                      interval(n).content.alphaL=(1-V(1))*(U(1)+...
                          Delta-1)/Delta+V(1);
                      interval(n).content.alphaR=V(2)*(1-U(4))*...
                          (1-U(1))/Delta+U(2);
                      interval(n).content.betaL=U(3)*(1-V(4))*...
                          (1-V(1))/Delta+V(3);
                      interval(n).content.betaR=(1-U(4))*(V(4)+...
                          Delta-1)/Delta+U(4);
                      interval(n).content.deltaL=(1-V(1))*U(5)/...
                          Delta+V(5)+(V(1)-1)*U(3)*V(6)/Delta;
                      interval(n).content.deltaR=(1-U(4))*V(6)/...
                          Delta+U(6)+(U(4)-1)*V(2)*U(5)/Delta;
                   end
               end
               for i=1:M
                   n=interval(i).leftchild;
                   if n>0 
                      x=interval(i).content.niuL;
                      y=interval(i).content.niuR;
                      z=interval(i).content.niu;  
                      u1=interval(n).content.alphaL;
                      u3=interval(n).content.betaL;
                      u5=interval(n).content.deltaL;
                      v4=interval(n+1).content.betaR;
                      v2=interval(n+1).content.alphaR;
                      v6=interval(n+1).content.deltaR;
                      interval(n).content.niuL=x;             
                      interval(n+1).content.niuR=y;
                      interval(n).content.niu=z;
                      interval(n+1).content.niu=z;
                      NIU=[1,v2;u3,1]\[y*(1-v4)-z*v6;x*(1-u1)-z*u5];
                      interval(n).content.niuR=NIU(1);
                      interval(n+1).content.niuL=NIU(2);
                   end
               end
               X=zeros(1,M);
               j=0;
               for i=1:M
                   if interval(i).leftchild==0 && interval(i).exist==1
                      X(j+1)=i;
                      j=j+1;
                   end
               end
               Num=j;
               Mat=[a-1;0;0;0];
               for i=1:Num
                   K=Ch(interval(X(i)).content.Interval);
                   A=P(interval,X(i),PhiL,PhiR,GL,GR);
                   Z1=zeros(10,1);
                   Z2=zeros(10,1);
                   Z3=zeros(10,1);
                   for j=1:10
                       Z1(j)=F(K(j));
                       Z2(j)=PhiL(K(j));
                       Z3(j)=PhiR(K(j));                      
                   end
                   sigma=A\Z1+A\Z2*interval(X(i)).content.niuL+...
                       A\Z3*interval(X(i)).content.niuR;
                   Mat1=zeros(4,10);
                   for j=1:10
                       Mat1(1,j)=K(j);
                       Mat1(2,j)=sigma(j);
                       Mat1(3,j)=X(i);
                   end
                   Mat=[Mat Mat1];
               end
               [~,idx]=sort(Mat(1,:));
               Mat=Mat(:,idx);
               JL=zeros(1,Num+1);
               JR=zeros(1,Num+1);
               for i=1:Num
                   n=Mat(3,10*i);
                   a1=interval(n).content.alphaL;
                   b1=interval(n).content.betaL;
                   d1=interval(n).content.deltaL;
                   n1=interval(n).content.niuL;
                   n2=interval(n).content.niuR;
                   JL(i+1)=JL(i)+d1+n1*a1+n2*b1;
                   m=Mat(3,10*(Num+1-i));
                   a2=interval(m).content.alphaR;
                   b2=interval(m).content.betaR;
                   d2=interval(m).content.deltaR;
                   m1=interval(m).content.niuL;
                   m2=interval(m).content.niuR;
                   JR(Num+1-i)=JR(Num+2-i)+d2+m1*a2+m2*b2;
               end
               for i=1:Num
                   Y=interval(Mat(3,10*i)).content.Interval;
                   K=Ch(Y);
                   Z1=zeros(10,1);
                   Z2=zeros(10,1);
                   for t=1:10
                       Z1(t)=GL(K(t));
                       Z2(t)=GR(K(t));
                   end
                   for j=1:10           
                       ur=GR(K(j))/s*(JL(i)+ ...
                       LIntegral(Y(1),K(j),Y(2),Z1.*(Mat(2,10*(i-1)+...
                       2:10*i+1))'))+GL(K(j))/s*(JR(i+1)+ ...
                       RIntegral(Y(1),K(j),Y(2),Z2.*(Mat(2,10*(i-1)+...
                       2:10*i+1))'));
                       Mat(4,10*(i-1)+1+j)=ur;
                   end                   
               end
               Er=Error(Mat,Mat2);
               Mat2=Mat;
               nn=size(Mat,2);
               N=M;
         end
         Ur=zeros(2,nn-1);
         Ur(1,:)=Mat(1,2:nn);
         Ur(2,:)=Mat(2,2:nn); 
         Leftendpoint=zeros(1,Num);
         for i=1:Num
             Leftendpoint(i)=interval(Mat(3,10*i)).content.Interval(1);
         end
end
