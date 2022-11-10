#### HELM Method ####
#### Generalized ####
## HELM solves for the following linear system:
## Y_{ik} * V_k(s) + I_0 = s * Conj(S) / V^hat_k(s)
## Conj(Y_{ik}) * V^hat_k(s) = s * S / V_k(s)
## With the reflection condition Conj(V_k(Conj(s))) = V^hat_k(s) at s = 1
## The above linear system can be solved order by order with respect to s
## For P-U Buses, an additional constraint is V_k(s) * V^hat_k(s) = V_ref ^2 + s(K - V_ref ^2)

library(latex2exp)
library(pracma)

### Parameters
Node_Num = 128                                                      ##Number of nodes, the first of which the reference bus
Line_Num = 128                                                      ##Number of lines
V_ref = 1                                                           ##Voltage of the reference bus
Z_ref = complex(real = 0, imaginary = 1) / 100                      ##Reference Impedence
Cond = diag(rep(1 / Z_ref, Line_Num))                               ##Conductance of the Lines
I_0 = rep(0, Node_Num)                                              ##Shunt current at each node
PQ_Bus = c(2:32, 34:64, 66:96, 98:128)                              ##Buses with P-Q constraints
PU_Bus = c(33, 65, 97)                                              ##Buses with P-U constraints

# Rand_Angle = rnorm(length(PQ_Bus), 0, 2*pi)
# Rand_Phase = complex(real = cos(Rand_Angle), imaginary = sin(Rand_Angle))
# Rand_radius = rnorm(length(PQ_Bus), 0, 35)
# S_PQ = -Rand_radius * Rand_Phase / 20                             ##Power Output at the P-Q buses
S_PQ = -c(rep(complex(real = 1, imaginary = .05), length(PQ_Bus))) / 20
P_PU = -c(rep(sum(Re(S_PQ)), length(PU_Bus))) / 4                   ##Real Power Output at the P-U buses
#U_PU = c(rep(1, length(PU_Bus)))                                   ##Voltage at the P-U buses
U_PU = c(1, 1, 1)

## Topology of the network
Topology = matrix(nrow = Line_Num, ncol = 2)
colnames(Topology) = c("Out", "In")

for(i in 1:(nrow(Topology)-1)){
  Topology[i, ] = c(i, i+1)
}
Topology[nrow(Topology), ] = c(nrow(Topology), 1)

# Topology in form of Node * Line
NL = matrix(0, ncol = Node_Num, nrow = Line_Num)
for(i in 1:Line_Num){
  NL[i, Topology[i, 1]] = 1
  NL[i, Topology[i, 2]] = -1
}

## Equivalent Admittance Matrix
Y_eq = t(NL) %*% Cond %*% NL

## Reordering of the Equivalent Admittance Matrix
# The order will be reference, P-Q, and P-U
Permutation = diag(rep(1, Line_Num))
Permutation = Permutation[,c(1, PQ_Bus, PU_Bus)]
Y_sorted = solve(Permutation) %*% Y_eq %*% Permutation

### Trivial Solution
# Solution when no current is fed to the system (open circuit solution)
V_0 = solve(Y_sorted[2:Node_Num, 2:Node_Num], -V_ref * Y_sorted[2:Node_Num, 1] - I_0[2:Node_Num])

### Coefficients for higher order terms
power_terms = 51                                          ##Number of Orders to be calculated
A = matrix(nrow = power_terms, ncol = Node_Num - 1)       ##Coefficients for V(s)
B = matrix(nrow = power_terms, ncol = Node_Num - 1)       ##Coefficients for 1/V^hat(s)
C = matrix(nrow = power_terms, ncol = Node_Num - 1)       ##Coefficients for V^hat(s)
D = matrix(nrow = power_terms, ncol = Node_Num - 1)       ##Coefficients for 1/V(s)
E = matrix(nrow = power_terms - 1, ncol = length(PU_Bus)) ##Coefficients for Q(s)

## Coefficients corresponding to the trivial solution
A[1,] = V_0
C[1,] = Conj(V_0)
B[1,] = 1 / C[1,]
D[1,] = 1 / A[1,]

## Matrix Used for Iteration Process
M = matrix(0, nrow = 4*(Node_Num - 1) + length(PU_Bus), ncol = 4*(Node_Num - 1) + length(PU_Bus))

# Equations for Y_{ik} * V_k(s) = s * Conj(S) / V^hat_k(s)
M[1:(Node_Num - 1), 1:(Node_Num - 1)] = Y_sorted[2:Node_Num, 2:Node_Num]
M[length(PQ_Bus) + 1:length(PU_Bus), 4*(Node_Num - 1) + 1:length(PU_Bus)] = complex(imaginary = 1) * diag(B[1, length(PQ_Bus) + 1:length(PU_Bus)])
# Equations for Conj(Y_{ik}) * V^hat_k(s) = s * S / V_k(s)
M[(Node_Num - 1) + 1:(Node_Num - 1), 2*(Node_Num - 1) + 1:(Node_Num - 1)] = Conj(Y_sorted[2:Node_Num, 2:Node_Num])
M[(Node_Num - 1) + length(PQ_Bus) + 1:length(PU_Bus), 4*(Node_Num - 1) + 1:length(PU_Bus)] = complex(imaginary = -1) * diag(D[1, length(PQ_Bus) + 1:length(PU_Bus)])
# Equations for V_k(s) / V_k(s) = 1
M[2*(Node_Num - 1) + 1:(Node_Num - 1), 1:(Node_Num - 1)] = diag(D[1,])
M[2*(Node_Num - 1) + 1:(Node_Num - 1), 3*(Node_Num - 1) + 1:(Node_Num - 1)] = diag(A[1,])
# Equations for V^hat_k(s) / V^hat_k(s) = 1
M[3*(Node_Num - 1) + 1:(Node_Num - 1), (Node_Num - 1) + 1:(Node_Num - 1)] = diag(C[1,])
M[3*(Node_Num - 1) + 1:(Node_Num - 1), 2*(Node_Num - 1) + 1:(Node_Num - 1)] = diag(B[1,])
# Equations for V_k(s) * V^hat_k(s) = 1 + s ((U_k) * Conj(U_k) - 1)
M[4*(Node_Num - 1) + 1:length(PU_Bus), length(PQ_Bus) + 1:length(PU_Bus)] = diag(C[1, length(PQ_Bus) + 1:length(PU_Bus)])
M[4*(Node_Num - 1) + 1:length(PU_Bus), 2*(Node_Num - 1) + length(PQ_Bus) + 1:length(PU_Bus)] = diag(C[1, length(PQ_Bus) + 1:length(PU_Bus)])

# Inverse of the Governing Matrix
start_time <- Sys.time()
N = qr.solve(M)
end_time <- Sys.time()
end_time - start_time

# start_time <- Sys.time()
# N = solve(M)
# end_time <- Sys.time()
# end_time - start_time

start_time <- Sys.time()
## Iteration Process
k = 1

rhs = c()
# Equations for Y_{ik} * V_k(s) = s * Conj(S) / V^hat_k(s)
rhs[1:length(PQ_Bus)] = Conj(S_PQ) * B[1, 1:length(PQ_Bus)]
rhs[length(PQ_Bus) + 1:length(PU_Bus)] = P_PU * B[1, length(PQ_Bus) + 1:length(PU_Bus)]
# Equations for Conj(Y_{ik}) * V^hat_k(s) = s * S / V_k(s)
rhs[(Node_Num - 1) + 1:length(PQ_Bus)] = S_PQ * D[1, 1:length(PQ_Bus)]
rhs[(Node_Num - 1) + length(PQ_Bus) + 1:length(PU_Bus)] = P_PU * D[1, length(PQ_Bus) + 1:length(PU_Bus)]
# Equations for V_k(s) / V_k(s) = 1
rhs[2*(Node_Num - 1) + 1:(Node_Num - 1)] = 0
# Equations for V^hat_k(s) / V^hat_k(s) = 1
rhs[3*(Node_Num - 1) + 1:(Node_Num - 1)] = 0
# Equations for V_k(s) * V^hat_k(s) = 1 + s ((U_k) * Conj(U_k) - 1)
rhs[4*(Node_Num - 1) + 1:length(PU_Bus)] = U_PU * Conj(U_PU) - (V_0 * Conj(V_0))[length(PQ_Bus) + 1:length(PU_Bus)]

# The solution to the linear system is the vector of coefficients of the embedded system
sol = N %*% rhs
A[2, ] = sol[1:(Node_Num - 1)]
B[2, ] = sol[(Node_Num - 1) + 1:(Node_Num - 1)]
C[2, ] = sol[2*(Node_Num - 1) + 1:(Node_Num - 1)]
D[2, ] = sol[3*(Node_Num - 1) + 1:(Node_Num - 1)]
E[1, ] = sol[4*(Node_Num - 1) + 1:length(PU_Bus)]

sum = rep(0, length(rhs))

tolerance = 10^(-10) # Tolerance Level for the Loop
for(k in 2:(power_terms-1)){
  rhs = c()
  # Equations for Y_{ik} * V_k(s) = s * Conj(S) / V^hat_k(s)
  rhs[1:length(PQ_Bus)] = Conj(S_PQ) * B[k, 1:length(PQ_Bus)]
  rhs[length(PQ_Bus) + 1:length(PU_Bus)] = P_PU * B[k, length(PQ_Bus) + 1:length(PU_Bus)]
  # Equations for Conj(Y_{ik}) * V^hat_k(s) = s * S / V_k(s)
  rhs[(Node_Num - 1) + 1:length(PQ_Bus)] = S_PQ * D[k, 1:length(PQ_Bus)]
  rhs[(Node_Num - 1) + length(PQ_Bus) + 1:length(PU_Bus)] = P_PU * D[k, length(PQ_Bus) + 1:length(PU_Bus)]
  # Equations for V_k(s) / V_k(s) = 1
  rhs[2*(Node_Num - 1) + 1:(Node_Num - 1)] = 0
  # Equations for V^hat_k(s) / V^hat_k(s) = 1
  rhs[3*(Node_Num - 1) + 1:(Node_Num - 1)] = 0
  # Equations for V_k(s) * V^hat_k(s) = 1 + s ((U_k) * Conj(U_k) - 1)
  rhs[4*(Node_Num - 1) + 1:length(PU_Bus)] = 0

  if(k >= 3){
    sum[length(PQ_Bus) + 1:length(PU_Bus)] = -complex(imaginary = 1) * apply(B[k:2, length(PQ_Bus) + 1:length(PU_Bus)] * E[1:(k-1), ], 2, sum)
    sum[(Node_Num - 1) + length(PQ_Bus) + 1:length(PU_Bus)] = complex(imaginary = 1) * apply(D[k:2, length(PQ_Bus) + 1:length(PU_Bus)] * E[1:(k-1), ], 2, sum)
    sum[2*(Node_Num - 1) + 1:(Node_Num - 1)] = -apply(A[k:2, ] * D[2:k, ], 2, sum)
    sum[3*(Node_Num - 1) + 1:(Node_Num - 1)] = -apply(B[k:2, ] * C[2:k, ], 2, sum)
    sum[4*(Node_Num - 1) + 1:length(PU_Bus)] = -apply(A[k:2, length(PQ_Bus) + 1:length(PU_Bus)] * C[2:k, length(PQ_Bus) + 1:length(PU_Bus)], 2, sum)
  }
  else{
    sum[length(PQ_Bus) + 1:length(PU_Bus)] = -complex(imaginary = 1) * B[2, length(PQ_Bus) + 1:length(PU_Bus)] * E[1, ]
    sum[(Node_Num - 1) + length(PQ_Bus) + 1:length(PU_Bus)] = complex(imaginary = 1) * D[2, length(PQ_Bus) + 1:length(PU_Bus)] * E[1, ]
    sum[2*(Node_Num - 1) + 1:(Node_Num - 1)] = -A[2, ] * D[2, ]
    sum[3*(Node_Num - 1) + 1:(Node_Num - 1)] = -B[2, ] * C[2, ]
    sum[4*(Node_Num - 1) + 1:length(PU_Bus)] = -A[2, length(PQ_Bus) + 1:length(PU_Bus)] * C[2, length(PQ_Bus) + 1:length(PU_Bus)]
  }
  
  # sum = rep(0, length(rhs))
  # for(n in 1:(k-1)){
  #   sum[length(PQ_Bus) + 1:length(PU_Bus)] = sum[length(PQ_Bus) + 1:length(PU_Bus)] - complex(imaginary = 1) * B[k+1-n, length(PQ_Bus) + 1:length(PU_Bus)] * E[n, ]
  #   sum[(Node_Num - 1) + length(PQ_Bus) + 1:length(PU_Bus)] = sum[(Node_Num - 1) + length(PQ_Bus) + 1:length(PU_Bus)] + complex(imaginary = 1) * D[k+1-n, length(PQ_Bus) + 1:length(PU_Bus)] * E[n, ]
  #   sum[2*(Node_Num - 1) + 1:(Node_Num - 1)] = sum[2*(Node_Num - 1) + 1:(Node_Num - 1)] - A[k+1-n, ] * D[n+1, ]
  #   sum[3*(Node_Num - 1) + 1:(Node_Num - 1)] = sum[3*(Node_Num - 1) + 1:(Node_Num - 1)] - B[k+1-n, ] * C[n+1, ]
  #   sum[4*(Node_Num - 1) + 1:length(PU_Bus)] = sum[4*(Node_Num - 1) + 1:length(PU_Bus)] - A[k+1-n, length(PQ_Bus) + 1:length(PU_Bus)] * C[n+1, length(PQ_Bus) + 1:length(PU_Bus)]
  # }

  rhs = rhs + sum

  sol = N %*% rhs
  A[k+1, ] = sol[1:(Node_Num - 1)]
  B[k+1, ] = sol[(Node_Num - 1) + 1:(Node_Num - 1)]
  C[k+1, ] = sol[2*(Node_Num - 1) + 1:(Node_Num - 1)]
  D[k+1, ] = sol[3*(Node_Num - 1) + 1:(Node_Num - 1)]
  E[k, ] = sol[4*(Node_Num - 1) + 1:length(PU_Bus)]

  if(max(abs(sol)) < tolerance * V_ref){
    break
  }
}
end_time <- Sys.time()
end_time - start_time

## Voltage at each node
V = apply(A, 2, sum, na.rm = TRUE)
V = Permutation %*% c(V_ref, V)
plot(abs(V), type = "l", xlab = TeX("Bus No."), ylab = TeX("Voltage (p.u.)"))

## Reactive Power at each P-U Bus
Q = apply(E, 2, sum, na.rm = TRUE)

# ## Sigma Plots
# sigma_R <- function(sigma_I){
#   sigma_I^2 - 1 / 4
# }
# Sigma = Conj(S_PQ) / diag(Cond)[PQ_Bus] / V_ref^2
# Sigma = c(Sigma, Conj(complex(real = P_PU, imaginary = Re(Q))) / diag(Cond)[PU_Bus] / V_ref^2)

# plot(seq(-1, 1, length.out = 100), sigma_R(seq(-1, 1, length.out = 100)), type = "l", xlab = TeX("$\\sigma_I$"), ylab = TeX("$\\sigma_R$"), asp = 1)
# lines(Im(Sigma), Re(Sigma), type = "p")
