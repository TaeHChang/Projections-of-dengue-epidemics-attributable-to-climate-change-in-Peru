##############################################################
########### Figure 3 AN projection visualization
############################################################

library(tidyr);library(dplyr);library(purrr) ;library(lubridate) ; library(readxl)


################################## 

# Changes in attributable numbers aggregated in 5-year intervals, compared to the baseline period.
# Results obtained from MCMC using 5,000 samples.
anrel_df <- read.csv("C:\\Users\\redwo\\OneDrive\\바탕 화면\\Peru Seminar/data/anrel_df.csv")

#################################################################################
### Figure 3. Time series of attributable numbers (AN). 
#### Time series graph showing changes in AN. 
#### Expressed as increases compared to the baseline period.

ggplot(anrel_df,  
       aes(x = Period, y = est, color = Scenario, group = Scenario)) +
  geom_line(position = position_dodge(width = 0.2), size = 2) + 
  geom_point(position = position_dodge(width = 0.2), size = 6) +  
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u),  # eCI 
                position = position_dodge(width = 0.2), width = 0.2) +  # Errorbar 
  scale_color_manual(values = c("#0080FF", "#FF9933", "#FF0000")) + 
  scale_y_continuous(limits = c(-35000, 100000)) +  
  labs(title = "Temperature-attributable cases by Scenario and Period",  
       x = "Period",  # x
       y = "Increase of cases",  # y
       color = "Scenario") +  
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black"))










