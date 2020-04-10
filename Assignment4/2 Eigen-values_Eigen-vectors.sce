clc;clear;close;
A=[]
disp("-------------------------------------------------------------------------")
printf("                             Enter a 3x3 matrix:\n\n")
for i=1:3
    for j=1:3
        printf("Enter element A(%d,%d):",i,j)
        A(i,j)=input("")
    end
end
disp("-------------------------------------------------------------------------")
disp(A,"The matrix is:")
l=poly(0,'lamda')
C=A-l*eye(3,3)
disp("-------------------------------------------------------------------------")
disp(C,"The characteristic matrix is:")
p=poly(A,'lambda')
disp(p,"The characteristic equation is:")
l=spec(A)
disp("-------------------------------------------------------------------------")
disp(l,"Eigen Values are:")
disp("-------------------------------------------------------------------------")
E=[]
C=[]
l=l'
for i=1:3
    B=A-l(i)*eye(3,3)
    printf("\n --> The characteristic matrix for lambda%d = %d is:\n",i,l(i))
    disp(B)
    R=rref(B)
    disp(R,"Row reduced form is: ")
    e=[]
    if(R(1,1)==0)then
        e(1,1)=1;e(2,1)=-R(2,1);e(3,1)=-R(3,1);
    elseif(R(2,2)==0)then
        e(1,1)=-R(1,2);e(2,1)=1;e(3,1)=-R(3,2);
    else
        e(1,1)=-R(1,3);e(2,1)=-R(2,3);e(3,1)=1;
    end
    c=e/norm(e)
    E=[E e]
    C=[C c]
end

disp("-------------------------------------------------------------------------")
disp("The Eigen Vectors are:")
disp(E)
disp("The Unit Eigen Vectors are:")
disp(C)
disp("-------------------------------------------------------------------------")
