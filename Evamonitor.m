function h=Evamonitor(str,n,phiL,phiR,gL,gR,F1)
         A=P(str,n,phiL,phiR,gL,gR);
         K=Ch(str(n).content.Interval);
         Z1=zeros(10,1);
         for i=1:10
             Z1(i)=F1(K(i));           
         end
         Theta=zeros(1,10);
         for i=1:10
             Theta(i)=(20-2*i+1)*pi/20;
         end
         CF=zeros(10,10);        
         for i=1:10
             for j=1:10
                 CF(i,j)=cos((i-1)*Theta(j))/5;
             end
         end
         CF(1,:)=CF(1,:)/2;
         H=CF*(A\Z1);
         h=abs(H(9))+abs(H(10)-H(8));
end