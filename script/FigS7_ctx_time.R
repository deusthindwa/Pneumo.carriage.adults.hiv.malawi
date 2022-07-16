A <-
left_join(
pcvpa.des %>%
  mutate(datex = year(date)) %>%
  group_by(datex) %>%
  tally() %>%
  rename("N" = n),

pcvpa.des %>%
  filter(ctx == "Yes") %>%
  mutate(datex = year(date)) %>%
  group_by(datex) %>%
  tally()) %>%
  
  mutate(prop = n/N) %>%

  ggplot() +
  geom_line(aes(x = datex, y = prop), size = 1) +
  theme_bw() +
  scale_y_continuous(limits = c(0.95, 1), labels = scales::percent_format(accuracy = 1)) + 
  labs(title = "", x = "Year of sampling", y = "Cotrimoxazole use") +
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 0, vjust = 0.5, hjust = 0.3), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold", size = 12)) +
  theme(legend.position = "none")

#combined plots
ggsave(here("output", "FigS7_ctx_time.png"),
       plot = (A),
       width = 6, height = 8, unit="in", dpi = 300)
