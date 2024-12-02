function B = wrpAngle(A)
A = A(:);
A = min(abs([A,A-360,A+360]),[],2);
A(A>180) = 360 - A(A>180);
B = A;
end