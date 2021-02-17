# creating interactive Plots for the output of homol.search by the nontarget package
# via plotly
library(plotly)
library(tibble)
library(ggplot2)

homo_tibble <- as_tibble(homol[[1]]) 

homo_tibble <- homo_tibble %>%
  mutate(hsid=homo_tibble$`HS IDs`) %>%
  mutate(mzsplit=homo_tibble$'m/z increment') %>%
  tidyr::separate_rows(c('hsid', 'mzsplit'), sep="/") %>%
  mutate(hsid=as.factor(hsid)) 

summary(homo_tibble)
homo_filter <- homo_tibble %>%
  filter(homo_tibble$"HS IDs" != "0")

homo_filter$mzsplit <- as.numeric(homo_filter$mzsplit)

homo_filter$col1 = NA
homo_filter$col1[homo_filter$mzsplit < 14] = " < 14"
homo_filter$col1[homo_filter$mzsplit < 15 & homo_filter$mzsplit > 14] = "> 14 x < 15"
homo_filter$col1[homo_filter$mzsplit < 18 & homo_filter$mzsplit > 15] = "> 15 x < 18"
homo_filter$col1[homo_filter$mzsplit < 21 & homo_filter$mzsplit > 18] = "> 18 x < 21"
homo_filter$col1[homo_filter$mzsplit < 30 & homo_filter$mzsplit > 21] = "> 21 x < 30"
homo_filter$col1[homo_filter$mzsplit < 32 & homo_filter$mzsplit > 30] = "> 30 x < 32"
homo_filter$col1[homo_filter$mzsplit < 90 & homo_filter$mzsplit > 32] = "> 32 x < 90"
homo_filter$col1[homo_filter$mzsplit < 200 & homo_filter$mzsplit > 90] = "> 90 x < 200"


homo_plot <- ggplot(data=homo_filter, aes(x = mz, y = RT, text = paste(
  "m/z increment ", homo_filter$`m/z increment`, "\n",
  "HS ID ", homo_filter$`HS IDs`, "\n",
  sep = ""))) +
  geom_point(size = 1, colour = "azure3")+
  geom_line(data=homo_filter, aes(x=mz, y=RT,group = hsid, color=col1, alpha = 0.4)) 


ggplotly(homo_plot)