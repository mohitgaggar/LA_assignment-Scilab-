clc;clear;close;
A=[]
disp("-------------------------------------------------------------------------")
printf("                          Enter a 3x3 matrix:\n\n")
for i=1:3
    for j=1:3
        printf("Enter element A(%d,%d):",i,j)
        A(i,j)=input("")
    end
end
disp(A,"The matrix is:")
x0=[]
disp("-------------------------------------------------------------------------")
printf("\nEnter the initial Eigen Vector:\n")
for i=1:3
    x0(i,1)=input("")
end
disp(x0,"Initial Eigen Vector")
a=max(x0)
disp(a,"Initial Largest Eigen Value")
disp("-------------------------------------------------------------------------")
v=A*x0
i=1
while abs(max(v)-a)>0.002 then
    printf("                         Iteration Number = %d\n",i)
    i=i+1
    a=max(v)
    disp(a,"Current Eigen Value:")
    x1=v/a
    disp(x1,"Current Eigen Vector:")
    v=A*x1
    disp("-------------------------------------------------------------------------")
end
format('v',5)
disp("-------------------------------------------------------------------------")
printf("                             Iteration Number = %d \n",i)
printf("             (Equal Eigen Vectors in iteration number %d and %d)\n",i-1,i)
disp("The largest Eigen Value:")
disp(max(v))
disp("The largest Eigen Vector:")
disp(v/a)
disp("-------------------------------------------------------------------------")
