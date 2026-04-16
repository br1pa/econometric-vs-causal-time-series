lars_structure   <- compute_bn_metrics(graph_lars, kb, df_full_graph)
print(lars_structure)
lar_structure    <- compute_bn_metrics(graph_lar, kb, df_full_graph)
print(lar_structure)
js_structure     <- compute_bn_metrics(graph_js, kb, df_full_graph)
print(js_structure)
simone_structure <- compute_bn_metrics(graph_simone, kb, df_full_graph)
print(simone_structure)
kb_structure     <- compute_bn_metrics(kb,kb, df_full_graph)
print(kb_structure)

shd_scores <- data.frame(
  Pair = c("LASSO vs LAR", "LASSO vs JS", "LASSO vs SIMONE", "LAR vs JS", "LAR vs SIMONE", "JS vs SIMONE"),
  SHD = c(
    shd(graph_lars, graph_lar),
    shd(graph_lars, graph_js),
    shd(graph_lars, graph_simone),
    shd(graph_lar, graph_js),
    shd(graph_lar, graph_simone),
    shd(graph_js, graph_simone)
  )
)

shd_scores <- shd_scores[order(shd_scores$SHD), ]

library(ggplot2)

# Plot the SHD scores
ggplot(shd_scores, aes(x = reorder(Pair, SHD), y = SHD, fill = Pair)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Pairwise SHD Score Comparison Among Learnt Structures",
       y = "SHD Score",
       x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))