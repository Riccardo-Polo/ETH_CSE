# My solution to exercise 1 of problem set 1
# In the first part we have computed analitycally the expressions 
# of the matrix F , the vector f and the constant c for the 
# unconstrained quadratic optimization problem. The expressions are:
# F = (calB)^T(calQ)(calB) + calR
# f = (calB)^T(calQ)(calA)x0
# c = (x0)^T Q x0+ (x0)^T (calA)^T (calQ) (calA) x0
# where
# calA = [A; A^2; ...; A^N]^T
# calB = [B; AB; ...; A^(N-1)B]^T
# calQ = diag(Q, Q, ..., Q,S) (T-1 times Q and 1 time S)
# calR = I_T times  R

# Then we have solved analytically the uncostrained problem getting the
# optimal input trajectory
# u* = -2 (F+F^T)^(-1) f^T

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_discrete_are

# The matrix A and B are given from the text and are
A = np.array([[0.88, 0.2, 0.49],
              [1.12, 0.94, -0.49],
              [0.48, -0.08, -0.05]])

B = np.array([[0.28, 0.36], [0.3, 0.27], [0.21, 0.33]])

# we define arbitrary Q,S,R matrices. For semplicity 
# we use diagonal matrices with positive entries on the diagonal.
Q = np.diag([1, 1, 1])
S = np.diag([1, 1, 1])
R = np.diag([1, 1])
x0 = np.array([1, 0, 0])

# solve the lqr problem with my closed formula and compare with a numerical solution obtained with scipy
N = 10
calA = np.zeros((N*3, 3))
calB = np.zeros((N*3, N*2))
for i in range(N):
    calA[i*3:(i+1)*3, :] = np.linalg.matrix_power(A, i+1)
    for j in range(i+1):
        calB[i*3:(i+1)*3, j*2:(j+1)*2] = np.linalg.matrix_power(A, i-j) @ B
calQ = np.zeros((N*3, N*3))
for i in range(N-1):
    calQ[i*3:(i+1)*3, i*3:(i+1)*3] = Q
calQ[(N-1)*3:N*3, (N-1)*3:N*3] = S
calR = np.zeros((N*2, N*2))
for i in range(N):
    calR[i*2:(i+1)*2, i*2:(i+1)*2] = R
F = calB.T @ calQ @ calB + calR
f = calB.T @ calQ @ calA @ x0
c = x0.T @ Q @ x0 + x0.T @ calA.T @ calQ @ calA @ x0
u_star = -2 * np.linalg.inv(F + F.T) @ f.T

# Now we solve the same problem with scipy
P = solve_discrete_are(A, B, Q, R)
K = np.linalg.inv(R + B.T @ P @ B) @ (B.T @ P @ A)
u_star_scipy = np.zeros((N*2,))
x = x0
for i in range(N):
    u_star_scipy[i*2:(i+1)*2] = -K @ x
    x = A @ x + B @ u_star_scipy[i*2:(i+1)*2]
# Now we compare the two solutions
print("Optimal input trajectory from my closed formula:")
print(u_star)
print("Optimal input trajectory from scipy:")
print(u_star_scipy)
# We can also plot the two trajectories to visually compare them #print one of the solutions as dotted lines and the other as solid lines
#print the dotted above the solid lines
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(u_star_scipy[0::2], label='u1 (scipy)', linestyle='-')
plt.plot(u_star[0::2], label='u1 (closed formula)', linestyle='--')
plt.title('Optimal input trajectory for u1')    
plt.xlabel('Time step')
plt.ylabel('Control input')
plt.legend()
plt.subplot(2, 1, 2)
plt.plot(u_star_scipy[1::2], label='u2 (scipy)', linestyle='-')
plt.plot(u_star[1::2], label='u2 (closed formula)', linestyle='--')
plt.title('Optimal input trajectory for u2')
plt.xlabel('Time step')
plt.ylabel('Control input')
plt.legend()
plt.tight_layout()
plt.show()


########################

import cvxpy as cp