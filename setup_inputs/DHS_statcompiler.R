library(haven)
library(rdhs)
library(survey)
library(data.table)
library(tidyverse)

options(survey.lonely.psu = "adjust")
if(!dir.exists('figures'))dir.create('figures')

#tagIds = c(10,32,33,36)
indicatorIDs <- c('ML_NETC_C_ITN', 'ML_AMLD_C_ACT', 'ML_PMAL_C_RDT', 'ML_PMAL_C_MSY')
dhs_df_raw <- dhs_data(indicatorIds = indicatorIDs, countryIds = c("GH"), breakdown = "subnational", surveyYearStart = 2000)
table(dhs_df_raw$SurveyId)
table(dhs_df_raw$CharacteristicLabel, dhs_df_raw$SurveyId)
#levels_to_keep <- grep("[..]", dhs_df_raw$CharacteristicLabel)
#dhs_df_raw <- dhs_df_raw[levels_to_keep,]
dhs_df_raw$NAME_1 <- gsub("[..]", "", dhs_df_raw$CharacteristicLabel)

table(dhs_df_raw$SurveyId, dhs_df_raw$Indicator)
dhs_df <- dhs_df_raw %>%
  filter(Indicator %in% c("Malaria prevalence according to microscopy", "Malaria prevalence according to RDT",
                          "Children who took any ACT", "Children under 5 who slept under an insecticide-treated net (ITN)")) %>%
  dplyr::select(SurveyYear, Indicator, Value, NAME_1) %>%
  mutate(SurveyYear = as.numeric(SurveyYear),
         Indicator = gsub("Malaria prevalence according to ", "", Indicator),
         Indicator = gsub("Children under 5 who slept under an insecticide-treated net ", "", Indicator),
         Indicator = gsub('[(]', '', gsub("[)]", "_U5", Indicator)),
         Indicator = gsub("Children who took any ACT", "ACT", Indicator)) #%>%
#pivot_wider(names_from = "Indicator", values_from = "Value")

pplot <- ggplot(data = subset(dhs_df)) +
  geom_point(aes(x = SurveyYear, y = Value, col = NAME_1)) +
  geom_line(aes(x = SurveyYear, y = Value, col = NAME_1)) +
  facet_wrap(~Indicator, scales = 'free', nrow = 2) +
  scale_x_continuous(lim = c(2005, 2020), breaks = seq(2005, 2020, 3), labels = seq(2005, 2020, 3)) +
  scale_y_continuous(lim = c(0, 100)) #+  scale_color_brewer(palette = 'Set2')

ggsave(file.path( 'figures','stat_compiler_GHA.png'), plot = pplot,   width = 12, height = 3.5, device = "png")

####EXAMPLE THIS NEEDS TO BE EDITED
dhs_df$ecozone <- NA
dhs_df$ecozone[dhs_df$NAME_1 %in% c('Ashanti','Central')] <- 'ecozone1'
dhs_df$ecozone[dhs_df$NAME_1 %in% c('Brong-Ahafo','Western')] <- 'ecozone2'
dhs_df$ecozone[dhs_df$NAME_1 %in% c('Greater Accra','Northern')] <- 'ecozone3'

pplot <- ggplot(data = subset(dhs_df)) +
  geom_point(aes(x = SurveyYear, y = Value, col = NAME_1)) +
  geom_line(aes(x = SurveyYear, y = Value, col = NAME_1)) +
  facet_wrap(ecozone ~ Indicator, scales = 'free', nrow = 2) +
  scale_x_continuous(lim = c(2005, 2020), breaks = seq(2005, 2020, 3), labels = seq(2005, 2020, 3)) +
  scale_y_continuous(lim = c(0, 100)) +
  scale_color_brewer(palette = 'Set2') +
  customTheme


ggsave(file.path( 'figures','stat_compiler_GHA_ecozone.png'), plot = pplot,   width = 12, height = 3.5, device = "png")
