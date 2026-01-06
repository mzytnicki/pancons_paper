library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_file_name  <- args[[1]]
col_id           <- as.integer(args[[2]])
output_file_name <- args[[3]]

input        <- read.table(input_file_name)
conservation <- input[, col_id]
start        <- input[, 2]
end          <- input[, 3]
size         <- end - start
conservation <- unlist(mapply(rep, conservation, size))

plot <- data.frame(conservation = conservation) |>
  ggplot2::ggplot(ggplot2::aes(x = conservation)) +
  ggplot2::geom_histogram() + 
  ggplot2::ylab("# base pairs") +
  ggplot2::theme_minimal()
ggsave(plot, file = output_file_name, height = 6, unit = "cm")
