library(dplyr)
library(ggplot2)
library(ggsignif)
library(data.table)
library(Rcmdr)

gastric <- CXXC5_Expression_Values[1:10,]
normal <- CXXC5_Expression_Values[10:20,]

normal_values <- c(normal$Average)
gastric_values <- c(gastric$Average)

n_avg <- mean(normal_values)
n_sd <- sd(normal_values)
n_var <- var(normal_values)
n_sem <- sqrt(n_var)/sqrt(length(normal_values))

n_lower <- n_avg - n_sem
n_upper <- n_avg + n_sem

g_avg <- mean(gastric_values)
g_sd <- sd(gastric_values)
g_var <- var(gastric_values)
g_sem <- sqrt(g_var)/sqrt(length(gastric_values))

g_lower <- g_avg - g_sem
g_upper <- g_avg + g_sem

lowers <- c(n_lower, g_lower)
uppers <- c(n_upper, g_upper)

expression_value <- c(n_avg, g_avg)

norm <- c("Normal")
norm <- rep(norm, length(normal_values))

gast <- c("Gastric Tumor")
gast <- rep(gast, length(gastric_values))

group <- c("Normal", "Gastric Tumor")

df = data.frame(group, expression_value)

df_errorbars = data.frame(group, expression_value, lowers, uppers)

p_value <- tibble(
  x = c("Normal", "Normal", "Gastric Tumor", "Gastric Tumor"),
  y = c(11.5, 11.5, 11.5, 11.5)
)

g <- ggplot(df, aes(x=group, y=expression_value)) + geom_bar(stat = "identity") +
  annotate("text", x=1.5, y=11.77, label="", size = 7) +
  coord_cartesian(ylim = c(0, 11.5)) + geom_errorbar(
    data=df_errorbars, mapping=aes(x=group, ymin=uppers, ymax=lowers,
                                   width = 0.3 )) + 
  labs(y="Mean CXXC5 Expression Level", x="Group")

g

u <- t.test(normal_values, gastric_values)
u

ggsave("cxxc5_gse79923.png")
