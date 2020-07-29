library(dplyr)
data <- read.csv('cdc_zika.csv')
data$value <- as.numeric(as.character(data$value))

a <- data %>% filter(report_date == date & data_field == 'cumulative_confirmed_local_cases')

b <- data %>% filter(location == loc & data_field == 'zika_reported_travel')

c <- data %>% top_n(1, value) %>%
  select(report_date, location)

d <- data %>% filter(location_type %in% types & data_field == 'zika_reported_travel') %>%
  select(report_date, location, data_field, value) %>%
  arrange(location)

e <- data %>% filter(location %in% c(loc, 'United_States-Florida') & data_field == 'zika_reported_travel')

f <- data %>% filter(location_type == 'county') %>%
  select(report_date, location, value) %>%
  top_n(1, value)

g <- data %>% filter(location == loc & report_date == date) %>%
  transmute(data_field, pct = 100*value/sum(value)) %>%
  arrange(desc(pct))

h <- data %>% filter(location_type == 'county') %>%
  group_by(data_field) %>%
  summarize(max_cases = max(value)) %>%
  arrange(data_field)

i <- data %>% filter(location_type == 'country') %>%
  group_by(location) %>%
  summarize(total = sum(value, na.rm=T)) %>%
  arrange(desc(total))

j <- data %>% filter(location_type == loc_type) %>%
  group_by(location) %>%
  summarize(total = sum(value, na.rm=T)) %>%
  transmute(location, pct = 100 * total / sum(total)) %>%
  arrange(pct)
