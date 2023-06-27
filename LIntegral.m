function I=LIntegral(u,x,v,Z) 
         I=0;
         CF=zeros(10,10); 
         Theta=zeros(1,10);
         for i=1:10
             Theta(i)=(20-2*i+1)*pi/20;
         end        
         for i=1:10
             for j=1:10
                 CF(i,j)=cos((i-1)*Theta(j))/5;
             end
         end
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
         L=(v-u)/2*L;
         J=L*CF*Z;
         for i=1:10
             I=I+J(i)*cos((i-1)*acos((2*x-u-v)/(v-u)));
         end
end