function Y=Ch1(X,K)       
         Y=zeros(K,1);
         for i=1:K
             Y(i)=(X(2)-X(1))/2*cos((2*K-2*i+1)*pi/(2*K))+(X(2)+X(1))/2;
         end
end