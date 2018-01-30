library(tidyverse)
library(flexsurv)

xx <- pgompertz(seq(0, 0.15, 0.001), shape = 1.77, rate = 1/0.03)
plot(200*seq(0, 0.15, 0.001), 1-xx, type = 'l')
yy <- pgompertz(seq(0, 0.15, 0.001), shape = 1.77, rate = (1/0.95)*(1/0.03))
lines(200*seq(0, 0.15, 0.001), 1-yy, col = "red")

zz <- pgompertz(seq(0, 0.15, 0.001), shape = 1.77, rate = (1/0.83)*(1/0.03))
lines(200*seq(0, 0.15, 0.001), 1-zz, col = "green")

abline(h = 0.2)
abline(v = 10)

abline(h = 0.4)
abline(v = 5)

yy <- pgompertz(seq(0, 0.15, 0.001), shape = 2.94, rate = 0.8*(1/0.02))
lines(seq(0, 0.15, 0.001), 1-yy, type = 'l', col = "red")

zz <- pgompertz(seq(0, 0.15, 0.001), shape = 2.94, rate = 1.2*(1/0.02))
lines(seq(0, 0.15, 0.001), 1-zz, type = 'l', col = "green")

rerun(100, {
  x <- rgompertz(n = 1000, shape = 2.94, rate = 1/0.02)
  lines(survfit(Surv(x, rep(1, 1000)) ~ 1), conf.int = F, col = "yellow")
})

lines(seq(0, 0.15, 0.001), 1-xx)

shape_pp <- 1/rnorm(100, mean = 0.02, sd = 0.008)
scale_pp <- rnorm(100, 2.94, 0.108)

xx <- pgompertz(seq(0, 0.15, 0.001), shape = 2.94, rate = 1/0.02)
plot(seq(0, 0.15, 0.001), 1-xx, type = 'n')

pwalk(list(x = shape_pp, y = scale_pp), function(x, y) {
  xx <- pgompertz(seq(0, 0.15, 0.005), shape = x, rate = y)
  lines(seq(0, 0.15, 0.005), 1-xx, type = 'l', col = "gray")
})

lines(seq(0, 0.15, 0.001), 1-xx, col = "black", lwd = 1.2)
#, xscale = 0.001

# PF state-----

xx <- pgompertz(seq(0, 0.15, 0.001), shape = 1.77, rate = 1/0.03)
plot(seq(0, 0.15, 0.001), 1-xx, type = 'l')

rerun(100, {
  x <- rgompertz(n = 1000, shape = 1.77, rate = 1/0.03)
  lines(survfit(Surv(x, rep(1, 1000)) ~ 1), conf.int = F, col = "gray")
})

x <- rgompertz(n = 1000, shape = 1.77, rate = 1/0.03)
plot(survfit(Surv(x, rep(1, 1000)) ~ 1), conf.int = F, col = "gray")

#ggsave("graph.pdf", ggplot(data = expand.grid(x = 0:200, y = 0:290), aes(x = x, y = y)) + geom_blank() + labs(x = "", y = "") + coord_equal() + theme_bw() + scale_x_continuous(breaks = seq(0, 200, 10), minor_breaks = 0:200) + scale_y_continuous(breaks = seq(0, 290, 10), minor_breaks = 0:290) + theme(panel.grid.major = element_line(colour = rgb(0.6, 0.6, 0.6)), panel.grid.minor = element_line(colour = "gray")), device = "pdf", width = 210, height = 297, units = "mm")

