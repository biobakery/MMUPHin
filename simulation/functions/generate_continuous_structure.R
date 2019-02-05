generate_continuous_structure <- function(nSample, nBatch) {
  batch <- sample.int(nBatch, size = nSample, replace = TRUE)
  tb_batch <- tibble::tibble(batch = batch) %>%
    dplyr::group_by(batch) %>%
    dplyr::summarise(n_batch = n())
  n_per_batch <- min(tb_batch$n_batch)
  continuous_structure <- runif(n = n_per_batch, min = -1, max = 1)
  batch <- purrr::map2(tb_batch$batch, tb_batch$n_batch, ~rep(.x, .y)) %>% unlist()
  score <- purrr::map(tb_batch$n_batch,
    ~ c(continuous_structure, runif(n = .x - n_per_batch, min = -1, max = 1))
  ) %>% unlist()
  data.frame(batch = as.factor(batch),
             score = score,
             row.names = paste0("Sample", 1:nSample))
}
