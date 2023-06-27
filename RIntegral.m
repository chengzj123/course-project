function I=RIntegral(u,x,v,Z) 
         I=0;
         CF=zeros(10,10); 
         Theta=zeros(1,10);
         for i=1:10
             Theta(i)=(20-2*i+1)*pi/20;
         end        
         for i=1:10
             for j=1:10
                 CF(i,j)=2*cos((i-1)*Theta(j))/10;
             end
         end
         CF(1,:)=CF(1,:)/2;
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
         R=(v-u)/2*R;
         J=R*CF*Z;
         for i=1:10
             I=I+J(i)*cos((i-1)*acos((2*x-u-v)/(v-u)));
         end
end
