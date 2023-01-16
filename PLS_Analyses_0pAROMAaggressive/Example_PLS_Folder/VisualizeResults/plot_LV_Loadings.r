library("dplyr")
library("ggplot2")
library("doBy")
library("BurStMisc")
library("wesanderson")
library("broom")
library("lme4")
library("scales")
library("Hmisc")



path = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

PLS_type="WBHC"

#LOAD RESULT SUMMARY FILES
PLS_Corrfile= paste("PLS", PLS_type, "LV_EffectCorrs.csv", sep="_")


PLS_corr_df<- read.csv(PLS_Corrfile, header=TRUE)


jpeg(paste("PLS", PLS_type, "LV1_CorrEffects.jpg",sep="_"), width = 600, height = 600)
PLS_corr_df %>%  ggplot(aes(Variable, LV1)) +
  geom_bar(stat = 'identity', position=position_dodge(), width = .5)+
  geom_errorbar(aes(Variable, ymin=LV1_LL, ymax=LV1_UL), width=0.05, position =position_dodge(0.7))+

  theme_minimal()+ theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),
                         legend.title.align = 0.5, legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=18), 
                         axis.text.y = element_text(size=18))+
  xlab("Age") + ylab("Correlation") + scale_y_continuous(position = "left")+
  ggtitle(label="LV1")+labs(subtitle="")
  # scale_fill_manual(name= "Effect Type",values= c(wes_palette(name="FantasticFox1")[2], wes_palette(name="FantasticFox1")[5]), breaks = c("Age", "Performance"), labels= c("Age", "Performance"))
dev.off()


jpeg(paste("PLS", PLS_type, "LV2_CorrEffects.jpg", sep="_"), width = 600, height = 600)
PLS_corr_df %>% ggplot(aes(Group, LV2)) +
  geom_bar(stat = 'identity', position=position_dodge(), width = .5)+
  geom_errorbar(aes(Group, ymin=LV2_LL, ymax=LV2_UL), width=0.05, position =position_dodge(0.7))+
  
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),
                         legend.title.align = 0.5, legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 18), axis.text.y = element_text(size = 15))+
  xlab("Group") + ylab("Correlations") + scale_y_continuous(position = "left" )+
  ggtitle(label="LV2 Correlations for Age")+labs(subtitle="")
# scale_fill_manual(name= "Effect Type",values= c(wes_palette(name="FantasticFox1")[2], wes_palette(name="FantasticFox1")[5]), breaks = c("Age", "Performance"), labels= c("Age", "Performance"))
dev.off()

