# inspired by https://github.com/rsquaredacademy/yahoofinancer/tree/master

fetch_prices <- function(ticker, interval,
                         start_date, end_date) {
  start <- as.numeric(as.POSIXct(start_date, tz = "UTC"))
  end <- as.numeric(as.POSIXct(end_date, tz = "UTC"))
  base_url <- "https://query2.finance.yahoo.com"
  path <- "v8/finance/chart/"
  url <- httr::modify_url(base_url,
    path = paste0(path, ticker)
  )

  query <- list(
    period1 = start,
    period2 = end,
    interval = interval
  )

  result <- httr::GET(url, query = query) |>
    httr::content("text", encoding = "UTF-8") |>
    jsonlite::fromJSON(simplifyVector = TRUE) |>
    (\(x) x$chart$result)()

  idx <- as.POSIXct(unlist(result$timestamp), origin = "1970-01-01", tz = "UTC")
  indicators <- result$indicators$quote[[1]]


  df <- data.frame(
    datetime = idx,
    close = unlist(indicators$close)
  )

  df
}

fetch_hourly_prices <- function(ticker = "AAPL",
                                start_date = "2020-01-01",
                                end_date = "2026-01-01") {
  fetch_prices(ticker, "1h", start_date, end_date)
}

fetch_daily_prices <- function(ticker = "AAPL",
                               start_date = "2020-01-01",
                               end_date = "2026-01-01") {
  fetch_prices(ticker, "1d", start_date, end_date)
}
