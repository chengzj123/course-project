function I=RIntegral1(u,x,v,Z,K) 
         I=0;
         CF=zeros(K,K); 
         Theta=zeros(1,K);
         for i=1:K
             Theta(i)=(2*K-2*i+1)*pi/(2*K);
         end        
         for i=1:K
             for j=1:K
                 CF(i,j)=2*cos((i-1)*Theta(j))/K;
             end
         end
         CF(1,:)=CF(1,:)/2;
         R=zeros(K,K);
         for i=3:K
             R(1,i)=(1/i-1/(i-2))/2;
             R(i-1,i-2)=-1/(2*(i-2));
             R(i-1,i)=1/(2*(i-2));
         end
         R(2,1)=-1;
         R(1,1)=1;
         R(1,2)=1/4;
         R(K,K-1)=-1/(2*(K-1));
         R=(v-u)/2*R;
         J=R*CF*Z;
         for i=1:K
             I=I+J(i)*cos((i-1)*acos((2*x-u-v)/(v-u)));
         end
end