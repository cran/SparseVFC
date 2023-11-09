## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(SparseVFC)
library(ggplot2)
library(dplyr)
library(tibble)

## -----------------------------------------------------------------------------
data(church)
X <- church$X
Y <- church$Y
CorrectIndex <- church$CorrectIndex

nX <- norm_vecs(X)
nY <- norm_vecs(Y)

## -----------------------------------------------------------------------------
set.seed(1614)
VecFld <- SparseVFC(nX, nY - nX, silent = FALSE)

## -----------------------------------------------------------------------------
vec <- expand.grid(x = seq(-1.2, 1.2, 0.2), y = seq(-1.2, 1.2, 0.2))
vec <- vec %>%
  rowwise() %>%
  mutate(v = list(predict(VecFld, c(x, y)))) %>%
  mutate(
    vx = v[1],
    vy = v[2]
  )

## -----------------------------------------------------------------------------
tibble(
  correct = 1:126 %in% CorrectIndex,
  VFC = 1:126 %in% VecFld$VFCIndex
) %>% table()

## -----------------------------------------------------------------------------
library(grid)
ggplot(vec, aes(x = x, y = y)) +
  geom_segment(aes(xend = x + vx, yend = y + vy),
    arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.25, alpha = 0.2
  ) +
  geom_segment(
    data = cbind(nX, nY - nX) %>% as.data.frame() %>% `colnames<-`(c("x", "y", "vx", "vy")),
    aes(xend = x + vx, yend = y + vy),
    arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.25
  ) +
  geom_segment(
    data = cbind(nX, nY - nX) %>% as.data.frame() %>% `colnames<-`(c("x", "y", "vx", "vy")) %>% slice(CorrectIndex),
    aes(xend = x + vx, yend = y + vy),
    arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.25, color = "red"
  )

