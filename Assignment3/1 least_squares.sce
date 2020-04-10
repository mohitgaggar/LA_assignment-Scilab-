clc;clear;close;
disp("-------------------------------------------------------------")
disp("              Projections by Least Squares")
disp("-------------------------------------------------------------")
n=input("Enter the number of linear quations of form C + Dt = b : ");
printf("\n")
A(1,1)=0
b(1,1)=0
max=0;
min=0;

//Creating the linear equations using input()
for(i=1:n)
    mprintf("Enter t%d: ",i);
    A(i,1)=input("")
    if(A(i,1)>max)
        max=A(i,1) ;
    end       
    if(A(i,1)<min)
        min=A(i,1);
    end
    mprintf("Enter b%d: ",i);
    b(i,1)=input("")
    printf("\n")
end
   
if(n==1)
    mprintf("\nThe best fit line is %d + %dt=b\n",A(1,1),b(1,1))
    return
end

a=A;
col_1=ones(n,1);
A=[col_1 A];
originalA=A
originalb=b;
disp(A,"A : ");
disp(b,"b: ")
clf();
scatter(A(:,2),b(:,1));
b=A'*b;
A=A'*A;
xhat=inv(A)*b
C=xhat(1,1)
D=xhat(2,1)
printf("\nThe value of C = %.2f\n",C)
printf("The value of D = %.2f\n",D)    
mprintf("The best fit line is %.2f + %.2f t = b\n",C,D)

disp("---------------------------------------------------")

disp(" Plotting the best fit line and original scatter plot : ")
x = min-3:1:max+3;
m = D;
c = C;
y = m * x + c;
plot2d(x, y,color("green"))

disp("---------------------------------------------------")

disp("VERIFICATION THAT ERROR IS MINIMUM")
p=originalA*xhat;
e=originalb-p;
printf("\nThe error vector e is : ")
disp(e)
printf("\nFor e to be minimum, a.e must be 0\n")
inner_product=sum(a.*e);
if(inner_product<0.00000000000001)
    inner_product=0;
end
mprintf("The dot/inner product a.e = %d\n",inner_product)

disp("-------------------END----------------------------------");
