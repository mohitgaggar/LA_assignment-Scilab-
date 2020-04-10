function m=makeMatrix(r, c)
// Creates a matrix from user input and returns it
  m = zeros(r, c)
  for i = 1:r
    for j = 1:c
      message = "Enter the element at [" + string(i) + "," + string(j) +"]"
      m(i,j) = input(message)
    end
  end
endfunction

function u = gaussianElimination(m)
// Performs gaussian elimination on given matrix, and returns upper triangular matrix
  u = m
  [r,c] = size(m)
  if c > r then
    c = r
  end
  for i = 1:c-1
    if u(i,i)~=0 then
      for j = i+1:r
        u(j,:) = u(j,:) - (u(j,i)/u(i,i))*u(i,:)
      end
    else
      flag = 0
      for j = i+1:r
        if u(j,i)~=0 then
          flag = 1
          for k = j:-1:i+1
            u([k,k-1],:) = u([k-1,k],:)
          end
        end
      end
      if flag~=0 then
        for j = i+1:r
          u(j,:) = u(j,:) - (u(j,i)/u(i,i))*u(i,:)
        end
      else
        disp("Breakdown of elimination occured.")
        break
      end
    end
  end
endfunction

function [p,l,u] = luDecomp(m)
// Performs LU decomposition of given matrix and returns the L and U matrices
  u = m
  [r,c] = size(m)
  l = eye(r, r)
  p = eye(r, r)
  for i = 1:c-1
    if u(i,i)~=0 then
      for j = i+1:r
        l(j,i) = u(j,i)/u(i,i)
        u(j,:) = u(j,:) - (u(j,i)/u(i,i))*u(i,:)
      end
    else
      flag = 0
      for j = i+1:r
        if m(j,i)~=0 then
          flag = 1
          for k = j:-1:i+1
            p([k,k-1],:) = p([k-1,k],:)
            u([k,k-1],:) = u([k-1,k],:)
          end
        end
      end
      if flag~=0 then
        for j = i+1:r
          l(j,i) = u(j,i)/u(i,i)
          u(j,:) = u(j,:) - (u(j,i)/u(i,i))*u(i,:)
        end
      else
        disp("Breakdown of elimination occured.")
        break
      end
    end
  end
endfunction

function n = invertMatrix(m)
// Returns inverse of matrix m
  [s,q] = size(m)
  if s!=q then
    disp("Inverse does not exist")
    abort
  end
  aug = gaussianElimination([m eye(s,s)])
  for i = 1:s
    if aug(i,i)==0 then
      disp("Inverse does not exist")
      abort
    end
    aug(i,:) = aug(i,:)/aug(i,i)
  end
  for i = s:-1:1
    for j = i-1:-1:1
      aug(j,:) = aug(j,:) - (aug(j,i)/aug(i,i))*aug(i,:)
    end
  end
  n = aug(:,s+1:2*s)

endfunction
© 2020 GitHub, Inc.
