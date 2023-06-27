function er=Error(A,B)
         n=size(A,2);
         m=size(B,2);
         e=0;
         f=0;
         for i=3:m
             for j=3:n
                 if A(1,j)>=B(1,i)
                    a=(B(1,i)-A(1,j-1))/(A(1,j)-A(1,j-1));
                    e=max(e,abs(a*A(4,j)+(1-a)*A(4,j-1)-B(4,i))); 
                    f=max(f,abs(a*A(4,j)+(1-a)*A(4,j-1)+B(4,i))); 
                    break
                 end
             end
         end
         er=e/f; 
end