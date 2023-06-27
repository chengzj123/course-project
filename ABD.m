function Z=ABD(str,n,F1,phiL,phiR,gL,gR)
         X=str(n).content.Interval;
         K=Ch(X);
         Z1=zeros(10,1);
         Z2=zeros(10,1);
         Z3=zeros(10,1);
         Z4=zeros(10,1);
         Z5=zeros(10,1);
         for i=1:10
             Z1(i)=F1(K(i));
             Z2(i)=phiL(K(i));
             Z3(i)=phiR(K(i));
             Z4(i)=gL(K(i));
             Z5(i)=gR(K(i));         
         end
         A=P(str,n,phiL,phiR,gL,gR);
         Z(1)=LIntegral(X(1),X(2),X(2),Z4.*(A\Z2));
         Z(2)=LIntegral(X(1),X(2),X(2),Z5.*(A\Z2));
         Z(3)=LIntegral(X(1),X(2),X(2),Z4.*(A\Z3));
         Z(4)=LIntegral(X(1),X(2),X(2),Z5.*(A\Z3));
         Z(5)=LIntegral(X(1),X(2),X(2),Z4.*(A\Z1));
         Z(6)=LIntegral(X(1),X(2),X(2),Z5.*(A\Z1));
end