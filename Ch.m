function K=Ch(X)       
         K=zeros(1,10);
         for i=1:10
             K(i)=(X(2)-X(1))/2*cos((20-2*i+1)*pi/20)+(X(2)+X(1))/2;
         end
end