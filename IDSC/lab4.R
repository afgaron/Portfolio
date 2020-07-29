# exercise 1
library(rvest)
library(dplyr)
nodes <- read_html('https://www.fdic.gov/bank/individual/failed/banklist.html')
tables <- html_nodes(nodes, 'table')
banklist <- html_table(tables[[1]])
count <- banklist %>% filter(City==city_st[1] & ST==city_st[2]) %>% 
  dplyr::count() %>% 
  as.numeric()
banks_df <- banklist %>% filter(ST %in% states) %>% 
  select(`Bank Name`)
banks <- banks_df$`Bank Name`

# exercise 2
nodes2 <- read_html('https://www.ohe.state.mn.us/dPg.cfm?pageID=792')
tables2 <- html_nodes(nodes2, 'table')
act_df <- html_table(tables2[[3]])[-1,1:3]
names(act_df) <- c('year', 'mn_composite', 'national_composite')
mn_better <- act_df[act_df$year==year,]$mn_composite > 
  act_df[act_df$year==year,]$national_composite

# exercise 3
deg_df <- read.csv('degrees.txt', header=F, col.names=c('degree', 'y1981', 'y2010'))
deg_df <- deg_df %>% mutate(pct_change = (y2010-y1981) / y1981 * 100) %>%
	filter(pct_change>pct) %>%
	arrange(desc(pct_change))

# exercise 4
states_df <- read.csv('us_states.csv')
states_df <- states_df %>% transmute(state, pop_density = population / area) %>% 
	filter(pop_density > min_density) %>% 
	arrange(pop_density)
