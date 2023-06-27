function M=P(str,n,phiL,phiR,gL,gR)        
         X=str(n).content.Interval;
         K=Ch(X);
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
         CB=5*CF';
         CF(1,:)=CF(1,:)/2;
         L=zeros(10,10);
         for i=3:10
             L(1,i)=(-1)^(i-1)*(1/i-1/(i-2))/2;
             L(i-1,i-2)=1/(2*(i-2));
             L(i-1,i)=-1/(2*(i-2));
         end
         L(2,1)=1;
         L(1,1)=1;
         L(1,2)=-1/4;
         L(10,9)=1/18;
         L=(X(2)-X(1))/2*L;
         R=zeros(10,10);
         for i=3:10
             R(1,i)=(1/i-1/(i-2))/2;
             R(i-1,i-2)=-1/(2*(i-2));
             R(i-1,i)=1/(2*(i-2));
         end
         R(2,1)=-1;
         R(1,1)=1;
         R(1,2)=1/4;
         R(10,9)=-1/18;
         R=(X(2)-X(1))/2*R;
         IL=CB*L*CF;
         IR=CB*R*CF;
         DPL=eye(10);
         DPR=eye(10);
         DGL=eye(10);     
         DGR=eye(10);
         for i=1:10
             DPR(i,i)=phiR(K(i));
             DGL(i,i)=gL(K(i));
             DGR(i,i)=gR(K(i));
             DPL(i,i)=phiL(K(i));            
         end
         M=eye(10)+DPL*IL*DGL+DPR*IR*DGR;      
end