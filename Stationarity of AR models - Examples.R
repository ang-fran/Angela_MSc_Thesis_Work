rm(list = ls())
T = 500

set.seed(0)

# AR (1) stat
y =numeric(T)
e = numeric(T)

y[1] = rnorm(1, 0, 1)
e[1] = rnorm(1, 0, 1)


for (t in 2:T) {
  e[t] = rnorm(1, 0, 1)
  y[t] = (0.7 * y[t - 1]) + e[t]
}

ts.plot(y)

# AR (1) stat
y = numeric(T)
e = numeric(T)

y[1] = rnorm(1, 0, 1)
e[1] = rnorm(1, 0, 1)


for (t in 2:T) {
  e[t] = rnorm(1, 0, 1)
  y[t] = (-1.01 * y[t - 1]) + e[t]
}

ts.plot(y)



# AR(2) - nonstat
T = 500

set.seed(0)
y =numeric(T)
e = numeric(T)

y[1] = rnorm(1, 0, 1)
y[2] = rnorm(1, 0, 1)
e[1] = rnorm(1, 0, 1)
e[2] = rnorm(1, 0, 1)

for (t in 3:T) {
  e[t] = rnorm(1, 0, 1)
  y[t] = (1.1 * y[t - 1]) + (0.8 * y[t - 2]) + e[t]
}

ts.plot(y)

# Stat AR(2)
T = 500

set.seed(0)
y =numeric(T)
e = numeric(T)

y[1] = rnorm(1, 0, 1)
y[2] = rnorm(1, 0, 1)
e[1] = rnorm(1, 0, 1)
e[2] = rnorm(1, 0, 1)

for (t in 3:T) {
  e[t] = rnorm(1, 0, 1)
  y[t] = (-0.2 * y[t - 1]) + (0.3 * y[t - 2]) + e[t]
}

ts.plot(y)

