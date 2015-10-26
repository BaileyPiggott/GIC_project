#libraries
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

# simulate data --------------------------
f <- 0.13333
n <- 1000
time <- 0:(n-1)
amp <- 2
signal <- amp*sin(2*pi*f*time) # no phase offset here
noise <- rnorm(n) # white noise
data <- signal + noise

plot(time, data, type='l') # plot of simulated data

# power spectrum of simulated data ---------------
periodogram <- abs(fft(data))^2
pgram <- periodogram[1:(n/2+1)] # spectra are symmetric, only need half to see what it looks like
freq <- seq(0, 0.5, by=1/n) # unit frequencies go from -0.5 to 0.5

# notice the frequency of the sinusoid is where we see a peak in the spectrum as well
plot(freq, pgram, type='l', log='y', main = 'no window')
abline(v=f, col='blue', lty=2)

# Add Hamming window ---------------------

# There are some windows on: https://en.wikipedia.org/wiki/Window_function
# try the Hamming window with alpha = 0.54 and beta = 0.46 : about 2/5 of the way down
alpha <- 0.54
beta <- 0.46
hamming <- alpha - beta*cos(2*pi*time/(n-1))

plot(time, data*hamming, type = 'l', main = 'windowed data') # data with window applied

specEst <- abs(fft(data*hamming))^2
pgram_est <- specEst[1:(n/2+1)]

plot(freq, pgram_est, type='l', log='y', main = 'with window')

# Geomag Data ------------------------------

# this is geomagnetic data from the Boulder geomagnetic observatory
# the sampling rate is 1 minute - this is likely more than we need, but we can worry about that later
load("boulder-geomag-Xdirection_1999-2012.RData")

x <- boulderX$value - mean(boulderX$value) # need to remove a mean from the data first

# only use a portion of these ~ 7 million points
N <- 100000
x <- x[1:N]

# a naive spectral estimate is the periodogram
# the (1/N) piece is actually a weighting of the data (equal weighting) - a "boxcar" taper (taper and window are the same thing in the literature)
P <- (1/N) * abs(fft(x))^2

# Trying plotting without the log='y' part - can't really see anything
plot(P, type='l', log='y')

# notice that the spectrum is symmetric, we only need half of it to see what's going on
P2 <- P[1:(N/2 + 1)] # we want the zeroth frequency and frequency at f=0.5

# To get actual frequencies on the x-axis:
freq <- seq(0, 0.5, by=1/N)

plot(freq, P2, type='l', log='y') 

# plot with ggplot ---------------
data <- data.frame(cbind(freq, P2)) %>% transmute(frequency = freq, periodogram = P2)

ggplot(data, aes(x = frequency, y = periodogram)) +
  geom_line() +
  scale_y_log10()


