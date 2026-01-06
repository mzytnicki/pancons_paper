args <- commandArgs(trailingOnly = TRUE)
inputFileNames <- head(args, -2)
outputFileName <- args[length(args) - 1]
ymax <- as.numeric(args[length(args)])

p <- purrr::map_dfr(inputFileNames, readr::read_tsv, col_names = c("score"), col_types = "n", n_max = 100000000, .id = "type") |>
  dplyr::ungroup() |>
  dplyr::mutate(type = factor(inputFileNames[as.numeric(type)])) |>
  dplyr::mutate(type = stringr::str_split_i(type, "/", -1)) |>
  dplyr::mutate(type = stringr::str_split_i(type, "_", -2)) |>
  ggplot2::ggplot(ggplot2::aes(x = score, colour = type)) +
  ggplot2::stat_ecdf(geom = "step", linetype = "dashed") +
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  ggplot2::ylim(0, ymax) +
  ggthemes::scale_color_colorblind() +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.title = ggplot2::element_blank())

ggplot2::ggsave(outputFileName, p)
