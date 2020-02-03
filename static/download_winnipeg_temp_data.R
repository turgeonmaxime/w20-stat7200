library(weathercan)
library(tidyverse)

# Download data--from Winnipeg The Forks
raw_data <- weather_dl(station_ids = 28051, 
                       start = "2010-01-01", 
                       end = "2019-12-31")

# Create wide dataset----
# Rows: Day of the year
# Cols: Year
clean_data <- raw_data %>% 
  group_by(year, date) %>% 
  summarise(avg_temp = mean(temp, na.rm = TRUE)) %>% # Average daily temp
  filter(month(date) != 2 | day(date) != 29) %>% # Remove leap years
  mutate(Date = paste(month(date), day(date), sep = "/")) %>% 
  select(-date) %>% 
  pivot_wider(names_from = year, names_prefix = "temp_",
              values_from = avg_temp) %>% 
  filter(!is.na(temp_2011), !is.na(temp_2012)) %>% # Remove rows with NAs
  select(-Date)

write_csv(clean_data,
          "static/winnipeg_temp.csv")
