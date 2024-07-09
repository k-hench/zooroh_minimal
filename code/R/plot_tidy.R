library(tidyverse)
library(here)
library(glue)
library(prismatic)
library(patchwork)
#library(ggrastr)
theme_set(theme_minimal() +
            theme(axis.line = element_line(linewidth = .25),
                  plot.subtitle = element_text(hjust = .5)))

data_f <- read_tsv(here("results/zooroh/inbreeding_by_HBD_class.tsv")) |> 
  pivot_longer(cols = -sample_id, 
               names_to = "HBDclass",
               values_to = "F",
               names_transform = as.numeric,
               names_prefix = "F_") |> 
  mutate(F = as.numeric(F))

data_roh <- read_tsv(here("results/zooroh/roh_segments.tsv.gz"))

data_local_3351 <- read_tsv(here("results/zooroh/local_hbd_prop/loacal_hbd_prop_ES3351.tsv.gz"))

hbd_classes <- sort(unique(data_f$HBDclass))

clrs <- c(rev(RColorBrewer::brewer.pal(ncol(data_local_3351)-2,"Set1")), "black") |> 
  set_names(nm = c(hbd_classes, "nonHBD"))

p1 <- data_f  |> 
  ggplot(aes(x = log2(HBDclass), 
             y = `F`, group = sample_id)) +
  geom_line(alpha = .3) +
  geom_point(alpha = .3) +
  labs(subtitle = "indv. proportion of genome in different HBD classes")

p2 <- data_f |> 
  group_by(sample_id, HBDclass) |> 
  summarise(cumulative_F = sum(F)) |> 
  mutate(F = cumulative_F - lag(cumulative_F,default = 0)) |> 
  ungroup() |> 
  mutate(HBDclass = factor(HBDclass, levels = rev(sort(unique(HBDclass))))) |> # filter(sample_id == "ES3351")
  ggplot(aes(x = sample_id, y = F, fill = HBDclass)) +
  geom_bar(stat = "identity",
           aes(color = after_scale(clr_darken(fill))),
           linewidth = .4) +
  scale_fill_manual(values = clrs,
                    labels = \(x){sprintf("%.0f",as.numeric(x))},
                    guide = guide_legend(nrow = 1)) +
  labs(subtitle = "Partitioning indiv. genomes in different HBD classes") +
  theme(axis.text.x = element_text(angle = 90))

chrms <- sort(unique(data_roh$chrom))

p3 <- data_roh |> 
  mutate(HBDclass = factor(HBDclass, levels = rev(sort(unique(HBDclass))))) |>
  filter(chrom == chrms[[1]]) |> # filter(id == "ES3351")
  ggplot(aes(y = id)) +
  geom_linerange(data = data_roh |> 
                   filter(!duplicated(id)),
                 aes(xmin = -Inf, xmax = Inf),
                 linewidth = .2) +
  geom_linerange(aes(xmin = start_pos, xmax = end_pos, color = HBDclass),
                 linewidth = 4) +
  labs(subtitle = "ROH Segments", y = "sample_id") +
  scale_color_manual(values = clrs, guide = "none") +
  scale_x_continuous(glue("Position on {chrms[[1]]}"),
                     labels = \(x){sprintf("%.0f cM", x)})

p4 <- data_local_3351 |> 
  pivot_longer(cols = -c(pos, sample_id)) |> 
    mutate(HBDclass = factor(str_remove(name, "^HBD_"),
           levels = c(sort(unique(data_f$HBDclass)), "nonHBD")))  |>  
  ggplot(aes(x = pos, y =  value, color = HBDclass)) +
  #  rasterize(geom_line(linewidth = .5), dpi = 300) +
  geom_line(linewidth = .5) +
  labs(subtitle = " By SNP for ES3351",
       y = "local HBD probabilities (=uncertainty)") +
  facet_grid(HBDclass ~ .) +
  scale_color_manual(values = clrs,
                     guide = "none")

layout <- "
ACD
BCD
BCD
"

pp <- p1 + p2 + p3 + p4 + 
  plot_layout(guides = "collect",
              design = layout) &
  theme(legend.position = "bottom")

ggsave(plot = pp, here("results/summary.pdf"), width = 14, height = 7, device = cairo_pdf)
