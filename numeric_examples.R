library(tidyverse)
library(cowplot)
library(ggthemes)

# seed variables
n <- 10
k <- 0.1
d <- seq(0.00001, 1, by = 0.001)

# generate data
df <- tibble(n, k, d)
consang <- df %>% 
  expand(n, k, d) %>% 
  mutate(h = (1 - d) / (1 - k*d)) %>% 
  mutate(q_dot = (h^2)/(n - (h^2) * (n - 1))) %>% 
  mutate(q = (1 / n) + ((n - 1)/n) * q_dot) %>% 
  mutate(cb_parent = ((q - ((1 - d)^2)*q)/(1 - ((1 - d)^2)*q))/2) %>% 
  mutate(cb_actdis = (q - ((1-d)^2)*q)/2) %>% 
  mutate(cb_actsed = ((q - ((1 - d)^2)*q)/(1 - (1 - d)*q))/2) %>% 
  mutate(cb_actdis_recdis = q/2) %>% 
  mutate(cb_actdis_recsed = (q - (1 - d)*q)/2) %>% 
  mutate(cb_actsed_recdis = ((q)/(1 - (1 - d)*q))/2) %>% 
  mutate(cb_actsed_recsed = ((q - (1 - d)*q)/(1 - (1 - d)*q))/2) %>% 
  mutate(cb_actdis_recdis_1 = 1/2) %>% 
  mutate(cb_actdis_recdis_0 = 0/2) %>% 
  mutate(cb_actdis_recsed_1 = (1 - (1 - d)*q)/2) %>% 
  mutate(cb_actdis_recsed_0 = -((1 - d)*q)/2) %>% 
  mutate(cb_actsed_recdis_1 = (1/(1 - (1 - d)*q))/2) %>%
  mutate(cb_actsed_recdis_0 = 0/2) %>%
  mutate(cb_actsed_recsed_1 = 1/2) %>%
  mutate(cb_actsed_recsed_0 = -((1 - d)*q)/(1 - (1 - d)*q)/2)

# build graphs
# generate theme
theme_cb <- theme_classic() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(margin = margin(t = 10, unit = "pt")),
    axis.text.y = element_text(margin = margin(r = 10)),
    axis.line = element_blank(),
    axis.ticks.length = unit(-5,"pt"),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

theme_cb_zoom <- theme_classic() +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(margin = margin(t = 10, unit = "pt")),
    axis.text.y = element_text(margin = margin(r = 10)),
    axis.line = element_blank(),
    axis.ticks.length = unit(-5,"pt"),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# dispersal and consanguinity
fudge_axis <- tibble(d = c(0,1), q = c(0.1,1))
f_consang <- consang %>% 
  ggplot(aes(x = d, y = q)) + geom_line() +  
  labs(x = expression(paste("Probability of dispersal (", italic("d"), ")")), 
       y = expression(paste("Neighborhood consanguinity (", italic("q"), ")"))) + 
  theme_cb_zoom + geom_rangeframe(data = fudge_axis) + scale_y_continuous(breaks = seq(0.1, 1, 0.1))
f_consang
save_plot(f_consang, path = "images", filename = "f_consang.pdf", base_height = 5, base_asp = 2)

# parent model
fudge_axis <- tibble(q = c(0.1,1), cb_parent = c(-2,2))
f_parent <- consang %>% 
  ggplot(aes(x = q, y = cb_parent)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
                              y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_parent

fudge_axis <- tibble(q = c(0.1,1), cb_parent = c(0.05,0.055))
f_parent_zoom <- ggplot(consang, aes(x = q, y = cb_parent)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_parent_zoom
save_plot(f_parent_zoom, path = "images", filename = "f_parent_zoom.pdf", base_height = 5, base_asp = 2)

# actor dispersing
fudge_axis <- tibble(q = c(0.1,1), cb_actdis = c(-2,2))
f_actdis <- ggplot(consang, aes(x = q, y = cb_actdis)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actdis

fudge_axis <- tibble(q = c(0.1,1), cb_actdis = c(0,0.05))
f_actdis_zoom <- ggplot(consang, aes(x = q, y = cb_actdis)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_actdis_zoom
save_plot(f_actdis_zoom, path = "images", filename = "f_actdis_zoom.pdf", base_height = 5, base_asp = 2)

# actor sedentary
fudge_axis <- tibble(q = c(0.1,1), cb_actsed = c(-2,2))
f_actsed <- ggplot(consang, aes(x = q, y = cb_actsed)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actsed

fudge_axis <- tibble(q = c(0.1,1), cb_actsed = c(0.05,0.058))
f_actsed_zoom <- ggplot(consang, aes(x = q, y = cb_actsed)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_actsed_zoom
save_plot(f_actsed_zoom, path = "images", filename = "f_actsed_zoom.pdf", base_height = 5, base_asp = 2)

# actor dispersing, recipient dispersing
fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recdis = c(-2,2))
f_actdis_recdis <- ggplot(consang, aes(x = q, y = cb_actdis_recdis)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actdis_recdis

fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recdis = c(0,0.5))
f_actdis_recdis_zoom <- ggplot(consang, aes(x = q, y = cb_actdis_recdis)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_actdis_recdis_zoom
save_plot(f_actdis_recdis_zoom, path = "images", filename = "f_actdis_recdis_zoom.pdf", base_height = 5, base_asp = 2)

# actor dispersing, recipient sedentary
fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recsed = c(-2,2))
f_actdis_recsed <- ggplot(consang, aes(x = q, y = cb_actdis_recsed)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actdis_recsed

fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recsed = c(0,0.05))
f_actdis_recsed_zoom <- ggplot(consang, aes(x = q, y = cb_actdis_recsed)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_actdis_recsed_zoom
save_plot(f_actdis_recsed_zoom, path = "images", filename = "f_actdis_recsed_zoom.pdf", base_height = 5, base_asp = 2)

# actor sedentary, recipient dispersing
fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recdis = c(-2,2))
f_actsed_recdis <- ggplot(consang, aes(x = q, y = cb_actsed_recdis)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  ylim(-2,2) + geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actsed_recdis

fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recdis = c(0,10))
f_actsed_recdis_zoom <- ggplot(consang, aes(x = q, y = cb_actsed_recdis)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis) + 
  ylim(0,10)
f_actsed_recdis_zoom
save_plot(f_actsed_recdis_zoom, path = "images", filename = "f_actsed_recdis_zoom.pdf", base_height = 5, base_asp = 2)

# actor sedentary, recipient sedentary
fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recsed = c(-2,2))
f_actsed_recsed <- ggplot(consang, aes(x = q, y = cb_actsed_recsed)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
        y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actsed_recsed

fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recsed = c(0.025,0.05))
f_actsed_recsed_zoom <- ggplot(consang, aes(x = q, y = cb_actsed_recsed)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_actsed_recsed_zoom
save_plot(f_actsed_recsed_zoom, path = "images", filename = "f_actsed_recsed_zoom.pdf", base_height = 5, base_asp = 2)

# actor dispersing, recipient dispersing, Q = 1
fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recdis_1 = c(-2,2))
f_actdis_recdis_1 <- ggplot(consang, aes(x = q, y = cb_actdis_recdis_1)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actdis_recdis_1

fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recdis_1 = c(0.45,0.55))
f_actdis_recdis_1_zoom <- ggplot(consang, aes(x = q, y = cb_actdis_recdis_1)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_actdis_recdis_1_zoom
save_plot(f_actdis_recdis_1_zoom, path = "images", filename = "f_actdis_recdis_1_zoom.pdf", base_height = 5, base_asp = 2)

# actor dispersing, recipient dispersing, Q = 0
fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recdis_0 = c(-2,2))
f_actdis_recdis_0 <- ggplot(consang, aes(x = q, y = cb_actdis_recdis_0)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actdis_recdis_0

fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recdis_0 = c(-0.05,0.05))
f_actdis_recdis_0_zoom <- ggplot(consang, aes(x = q, y = cb_actdis_recdis_0)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_actdis_recdis_0_zoom
save_plot(f_actdis_recdis_0_zoom, path = "images", filename = "f_actdis_recdis_0_zoom.pdf", base_height = 5, base_asp = 2)

# actor dispersing, recipient sedentary, Q = 1
fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recsed_1 = c(-2,2))
f_actdis_recsed_1 <- ggplot(consang, aes(x = q, y = cb_actdis_recsed_1)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actdis_recsed_1

fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recsed_1 = c(0,0.5))
f_actdis_recsed_1_zoom <- ggplot(consang, aes(x = q, y = cb_actdis_recsed_1)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_actdis_recsed_1_zoom
save_plot(f_actdis_recsed_1_zoom, path = "images", filename = "f_actdis_recsed_1_zoom.pdf", base_height = 5, base_asp = 2)

# actor dispersing, recipient sedentary, Q = 0
fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recsed_0 = c(-2,2))
f_actdis_recsed_0 <- ggplot(consang, aes(x = q, y = cb_actdis_recsed_0)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actdis_recsed_0

fudge_axis <- tibble(q = c(0.1,1), cb_actdis_recsed_0 = c(0,-0.5))
f_actdis_recsed_0_zoom <- ggplot(consang, aes(x = q, y = cb_actdis_recsed_0)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_actdis_recsed_0_zoom
save_plot(f_actdis_recsed_0_zoom, path = "images", filename = "f_actdis_recsed_0_zoom.pdf", base_height = 5, base_asp = 2)

# actor sedentary, recipient dispersing, Q = 1
fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recdis_1 = c(-2,2))
f_actsed_recdis_1 <- ggplot(consang, aes(x = q, y = cb_actsed_recdis_1)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  ylim(-2,2) + geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actsed_recdis_1

fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recdis_1 = c(0,10))
f_actsed_recdis_1_zoom <- ggplot(consang, aes(x = q, y = cb_actsed_recdis_1)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis) + 
  ylim(0,10)
f_actsed_recdis_1_zoom
save_plot(f_actsed_recdis_1_zoom, path = "images", filename = "f_actsed_recdis_1_zoom.pdf", base_height = 5, base_asp = 2)

# actor sedentary, recipient dispersing, Q = 0
fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recdis_0 = c(-2,2))
f_actsed_recdis_0 <- ggplot(consang, aes(x = q, y = cb_actsed_recdis_0)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actsed_recdis_0

fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recdis_0 = c(-0.05,0.05))
f_actsed_recdis_0_zoom <- ggplot(consang, aes(x = q, y = cb_actsed_recdis_0)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_actsed_recdis_0_zoom
save_plot(f_actsed_recdis_0_zoom, path = "images", filename = "f_actsed_recdis_0_zoom.pdf", base_height = 5, base_asp = 2)

# actor sedentary, recipient sedentary, Q = 1
fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recsed_1 = c(-2,2))
f_actsed_recsed_1 <- ggplot(consang, aes(x = q, y = cb_actsed_recsed_1)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25)
f_actsed_recsed_1

fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recsed_1 = c(0.45,0.55))
f_actsed_recsed_1_zoom <- ggplot(consang, aes(x = q, y = cb_actsed_recsed_1)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis)
f_actsed_recsed_1_zoom
save_plot(f_actsed_recsed_1_zoom, path = "images", filename = "f_actsed_recsed_1_zoom.pdf", base_height = 5, base_asp = 2)

# actor sedentary, recipient sedentary, Q = 0
fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recsed_0 = c(-2,2))
f_actsed_recsed_0 <- ggplot(consang, aes(x = q, y = cb_actsed_recsed_0)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb + geom_rangeframe(data = fudge_axis) + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + 
  geom_ribbon(aes(ymin = -1, ymax = 1), alpha=0.25) + ylim(-2,2)
f_actsed_recsed_0

fudge_axis <- tibble(q = c(0.1,1), cb_actsed_recsed_0 = c(0,-10))
f_actsed_recsed_0_zoom <- ggplot(consang, aes(x = q, y = cb_actsed_recsed_0)) + geom_line() +  
  labs(x = expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
       y = expression(paste("Evolutionarily stable ", italic("C/B")))) + 
  theme_cb_zoom + scale_x_continuous(breaks = seq(0.1, 1, 0.1)) + geom_rangeframe(data = fudge_axis) + 
  ylim(-10,0)
f_actsed_recsed_0_zoom
save_plot(f_actsed_recsed_0_zoom, path = "images", filename = "f_actsed_recsed_0_zoom.pdf", base_height = 5, base_asp = 2)

# generate combined plot of actor-only actual cost-benefit ratios (unscaled)
zoom_act <- plot_grid(f_actdis_zoom, f_actsed_zoom)
zoom_act

# generate combined plot of actual cost-benefit ratios on the same scale
scaled <- plot_grid(f_parent + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),  
          f_actdis + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),  
          f_actsed + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
          f_actdis_recdis + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
          f_actdis_recsed + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
          f_actsed_recdis + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
          f_actsed_recsed + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
          f_actdis_recdis_1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
          f_actdis_recdis_0 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
          f_actdis_recsed_1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
          f_actdis_recsed_0 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
          f_actsed_recdis_1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),  
          f_actsed_recdis_0 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
          f_actsed_recsed_1 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
          f_actsed_recsed_0 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
          ncol = 3, labels = "auto", label_x = 0.2, hjust = 0, label_size = 18, align = "hv")

y_axis <- ggdraw() + draw_label(expression(paste("Evolutionarily stable ", italic("C/B"))),
                                fontface = 'bold', size = 20, hjust = 0.5, angle = 90)
x_axis <- ggdraw() + draw_label(expression(paste("Neighborhood consanguinity (", italic("q"), ")")), 
                                fontface = 'bold', size = 20, hjust = 0.5)
cb_scaled <- plot_grid(y_axis, scaled, NULL, x_axis, nrow = 2, ncol = 2, rel_widths = c(0.1, 1), rel_heights = c(1, 0.1))
save_plot(cb_scaled, path = "images", filename = "cb_scaled.pdf", base_width = 8.5, base_height = 8)