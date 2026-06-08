# inspired by https://github.com/rsquaredacademy/yahoofinancer/tree/master

fetch_prices <- function(ticker, interval,
                         start_date, end_date) {
  # convert the dates to the required format
  start <- as.numeric(as.POSIXct(start_date, tz = "UTC"))
  end <- as.numeric(as.POSIXct(end_date, tz = "UTC"))

  # specify uri
  base_url <- "https://query2.finance.yahoo.com"
  path <- "v8/finance/chart/"
  url <- httr::modify_url(base_url,
    path = paste0(path, ticker)
  )
  
  # compose a query
  query <- list(
    period1 = start,
    period2 = end,
    interval = interval
  )

  # a file-based cache
  filename <- paste0("./data/", ticker, "-", start_date, "-", end_date, ".csv")
  if (file.exists(filename)) {
    df <- read.csv(filename)
    df$datetime <- as.POSIXct(df$datetime, tz = "UTC")
    return(df)
  } 
  
  # if not cached, fetch
  result <- httr::GET(url, query = query) |>
    httr::content("text", encoding = "UTF-8") |>
    jsonlite::fromJSON(simplifyVector = TRUE) |>
    (\(x) x$chart$result)()

  # convert 
  idx <- as.POSIXct(unlist(result$timestamp), 
    origin = "1970-01-01", tz = "UTC")
  indicators <- result$indicators$quote[[1]]
  
  df <- data.frame(
    datetime = idx,
    close = unlist(indicators$close)
  )

  write.csv(df, file = filename)
  df
}

fetch_hourly_prices <- function(ticker = "AAPL") {
  # only the last 730 days are available, so returns all of them
  today <- Sys.Date()
  start_date <- as.character(today - 720)
  end_date <- as.character(today)
  fetch_prices(ticker, "1h", start_date, end_date)
}

fetch_daily_prices <- function(ticker = "AAPL",
                               start_date = "2020-01-01",
                               end_date = "2026-01-01") {
  fetch_prices(ticker, "1d", start_date, end_date)
}
