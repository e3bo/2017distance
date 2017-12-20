#!/usr/bin/Rscript

library(stringr)
library(magrittr)

lines <- readLines("noise-map.tex")

move_label <- function(lns, label, dx = 0, dy = 0){
    is_lab_ln <- str_detect(lns, label)
    lab_lns <- lines[is_lab_ln]
    coord_pattern <- "(\\d+(\\.\\d+)?),\\s*(\\d+(\\.\\d+))"
    coords <- str_extract(lab_lns, coord_pattern)
    getc <- function(coord, pos) as.numeric(str_split(coord, ",")[[1]][pos])
    xcoords <- sapply(coords, getc, pos = 1)
    ycoords <- sapply(coords, getc, pos = 2)
    new_xcoords <- xcoords + dx
    new_ycoords <- ycoords + dy
    new_coords <- str_c(new_xcoords, new_ycoords, sep=",")
    new_lab_lns <- str_replace(lab_lns, coords, new_coords)
    new_lines <- lns
    new_lines[is_lab_ln] <- new_lab_lns
    new_lines
}

dx <- -10
new_lines <- lines %>% move_label("d_\\{11\\}", dx, dx) %>%
  move_label("d_\\{22\\}", -dx, dx) %>%
  move_label("d_\\{12\\}", dx, 0)

dx <- -5
new_lines <- new_lines %>% move_label("\\{-20\\};", dx) %>%
  move_label("\\{0\\};", dx) %>% move_label("\\{20\\};", dx) %>%
  move_label("\\{-10\\};", dx) %>% move_label("\\{10\\};", dx) %>%
  move_label("\\{-1\\};", dx) %>% move_label("\\{1\\};", dx)

dx <- -7.5
new_lines <- new_lines %>% move_label("\\{0.00\\};", -5, -5) %>%
  move_label("\\{0.05\\};", dx) %>% move_label("\\{0.10\\};", dx)

dx <- 7.5
new_lines <- new_lines %>% move_label("\\{0.000\\};", dx) %>%
  move_label("\\{0.005\\};", dx) %>% move_label("\\{0.010\\};", dx)

dx <- -5
new_lines <- new_lines %>% move_label("\\{0.0\\};", dx) %>%
  move_label("\\{0.2\\};", dx) %>% move_label("\\{0.4\\};", dx) %>%
  move_label("\\{0.02\\};", dx) %>% move_label("\\{0.04\\};", dx) %>%
  move_label("\\{0.02\\};", dx)

new_lines <- new_lines %>% move_label("211.38\\) \\{0.00\\};", -7.50, 0)

writeLines(new_lines, con="noise-map-adjusted.tex")
