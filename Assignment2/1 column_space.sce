clc;clear;close;

disp("---------------------------------------------------------------------")

    mprintf("\n")
    disp("           Column Space of a given matrix : ")
    mprintf("\n\nEnter the number of rows : \n\n")
    m=input("m= ")
    mprintf("\nEnter the number of columns : \n")
    n=input("n= ")
    
    mprintf("Enter values for matrix A :\n\n")
    for i=1:m
       for j=1:n
           mprintf("Enter value for A(%d,%d) = ",i,j)
           a(i,j)=input("")
       end
    end
    
    disp(a,"A:")
    
    a_temp=a
    M=m
    N=n
    
    k=1
    j_start=1;
    i=1;
    A=a_temp
    pivots_i=1;
    pivots(1,1)=0;
    while(i<=N)
        flag=0
        temp=0;
        for(j=j_start:M)
            if(A(j,i)<>0)
                flag=j;
                temp=i;
                break;
             end
        end
        if(flag<>0 )
            while(flag<>k)
                    A([flag,flag-1],:)=A([flag-1,flag],:);
                    flag=flag-1;
            end
            k=k+1;
            j_start=k;
            pivots(pivots_i,1)=flag;
            pivots(pivots_i,2)=temp;
            pivots_i=pivots_i+1;
            
            for(t=flag+1:M)
                 A(t,:)=A(t,:)-(A(t,temp)/A(flag,temp))*A(flag,:);
                 for(zx=1:N)
                     if(A(t,zx)<0.0000000000000001 & A(t,zx)>0)
                         A(t,zx)=0
                     end
                 end 
            end
        end
        i=i+1
     end
     Pivots=pivots
     Pivots_i=pivots_i
     a=A
     disp(A,"U: ")    // U
     mprintf("\n\nThe pivot columns are : ")
     for(v=1:(pivots_i-1))
         printf(" Col %d, ",pivots(v,2));
     end

     mprintf("\n\nColumn space of A has the Basis : {( ")
     for(v=1:(pivots_i-1))
         x=pivots(v,2);
         //disp(a_temp(:,x))
         for(w=1:M)
             mprintf("%.2f ,",a_temp(w,x));
         end
         if(v<>pivots_i-1)
            mprintf(" ),( ");
         end
     end       
     mprintf(")}")
     mprintf("\n\nRank of the matrix = %d\n",pivots_i-1)
     mprintf("\nDimension of column space = %d\n ",pivots_i-1);
disp("--------------------------------------------------------------------")
