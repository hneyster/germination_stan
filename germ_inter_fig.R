library(bayesplot)
library(bayestestR)
library(ggplot2)
library(cowplot)
ci_fxn <- function(x) {
  y <- as.data.frame(ci(x,method="HDI", ci=0.50)) #gives 89% CrI (using HDI)
  y[1] <- median(x)
  return(y)
}
xlabel <- 'Stratification (days)'
load("/home/harold/Dropbox/gitfiles/germ/mod_time_pois.rdata") #this is my rstanarm  stanfit object.
load("/home/harold/Dropbox/gitfiles/germ/mod_rate.rdata") #this nd the next models are  rstanarm stanfit objects
load("/home/harold/Dropbox/gitfiles/germ/mod_gr.rdata")

m <- mod_time_pois %>% as.matrix() %>% .[,1:16] %>% as.data.frame()


eu_s <- (m$`(Intercept)` + m$temp3 + m$strat +  m$`strat:temp3`) %>% ci_fxn
  us_s <- (m$`(Intercept)`+ m$origin + m$temp3 + m$strat + m$`origin:temp3` + m$`origin:strat` + m$`strat:temp3`+ m$`origin:strat:temp3`) %>% ci_fxn
  
eu <- (m$`(Intercept)` + m$temp3 %>% ci_fxn)
us <- (m$`(Intercept)`+ m$origin + m$temp3 + m$`origin:temp3`)  %>% ci_fxn
  
pp <- rbind(eu, us, eu_s, us_s)
pp <- cbind(pp, 'origin' = factor(c('native','introduced','native','introduced'), levels = c('native','introduced')))
pp <- cbind(pp, 'strat' = c("60","60","30","30"))
time_inter <- ggplot(pp, aes(x = strat, y = CI, ymin = CI_low, ymax = CI_high, color = origin, group = origin)) + geom_path() +
  geom_pointrange( position=position_jitter(width=0.08,height=0),size=1, alpha = 1) + scale_x_discrete( limits =c('60','30')) + 
  ylab("Germination timing at 32°C (log days)")+ xlab("") + xaxis_text(FALSE)+ legend_none()

### GROWTH ### 
m <- mod_gr %>% as.matrix() %>% .[,1:16] %>% as.data.frame()

eu_s <- (m$`(Intercept)` + m$temp2 + m$strat +  m$`strat:temp2`) %>% ci_fxn
us_s <- (m$`(Intercept)`+ m$origin + m$temp2 + m$strat + m$`origin:temp2` + m$`origin:strat` + m$`strat:temp2`+ m$`origin:strat:temp2`) %>% ci_fxn

eu <- (m$`(Intercept)` + m$temp2) %>% ci_fxn
us <-( m$`(Intercept)`+ m$origin + m$temp2 + m$`origin:temp2`)  %>% ci_fxn

pp <- rbind(eu, us, eu_s, us_s)
pp <- cbind(pp, 'origin' = factor(c('native','introduced','native','introduced'), levels = c('native','introduced')))
pp <- cbind(pp, 'strat' = c("60","60","30","30"))
gr_inter <-  ggplot(pp, aes(x = strat, y = CI, ymin = CI_low, ymax = CI_high, color = origin, group = origin)) + geom_path() +
  geom_pointrange( position=position_jitter(width=0.08,height=0),size=1, alpha = 1) + scale_x_discrete( limits =c('60','30')) +
  ylab("Growth rate at 27.3°C (cm/day)") + xlab(xlabel) + theme(legend.position = 'top') 

legend <- get_legend(
  # create some space to the left of the legend
  gr_inter + theme(legend.box.margin = margin(0, 0, 0, 20))
)
p0 <- plot_grid(time_inter, gr_inter +theme(legend.position = 'none'), ncol =1, rel_heights = c(1,1.1))
p1 <- plot_grid(legend, p0,rel_heights = c(.1,1), nrow = 2)
ggsave("germ_inter_fig.pdf",p1,width = 5,device ='pdf', height=7, units="in")

