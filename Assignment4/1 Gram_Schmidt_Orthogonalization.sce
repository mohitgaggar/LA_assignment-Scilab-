clc;clear;close;
A=[]
disp("-------------------------------------------------------------------------")
printf("\n           Enter 3 independent vectors to form of 3X3 matrix:\n")
for j=1:3
    printf("\n Enter vector v%d:\n",j);
    for i=1:3
       A(i,j)=input("");
    end
end
v1=A(:,1)
v2=A(:,2)
v3=A(:,3)
printf("The three independent vectors are:\n")
disp(v1,"v1=")
disp(v2,"v2=")
disp(v3,"v3=")
disp("-------------------------------------------------------------------------")
disp("Corresponding Orthonormal Vectors obtained using Gram Schimdt Process :")
q1=v1/norm(v1)
e=v2-(q1'*v2)*q1
q2=e/norm(e)
E=v3-((q1'*v3)*q1+(q2'*v3)*q2)
q3=E/norm(E)
disp(q1,"q1=")
disp(q2,"q2=")
disp(q3,"q3=")
disp("-------------------------------------------------------------------------")
