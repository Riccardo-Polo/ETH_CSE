clc
clear all

% write down A,B,C,D
A = [4/3 2/3;1 0];
B = [1 ; 0];
C = [-2/3 1];

Q = C.' * C + 0.001 * eye(2);
R = 0.001;
P_N = Q;

N = 5;
%% Point 1

x_0 = [10 ; 10];

[r,c] = size(P_N);
P = zeros(r,c*N);
P(1:r+1,1:c+1) = P_N;

[ru , cu] = size(B);
u_opt = zeros(ru,cu*(N-1));

for i=1:N
    P_next = P(1:r+1,c+2 : c+2+c);
    P_now = A.' * P_next * A + Q - A.' * P_next * B * inv(B.' * P_next *B + R)*B.' * P_next * A;

    u_opt(1:r,(i-1))
end
