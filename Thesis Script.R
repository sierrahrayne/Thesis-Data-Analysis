column_sum <- sum(density_df$density)
column_sum

#Figure 1
abundance_df$totalabun <- rowSums( abundance_df[,2:9] )
density_df <- abundance_df %>%
  mutate(density = (totalabun*2)) %>%
  mutate(MajorHabitat = ifelse(MajorHabitat == "Seagrass", "Seagrass Fringe", MajorHabitat))

kruskal.test(density ~ MajorHabitat, data = density_df)
dunn.test(density_df$density, g = density_df$MajorHabitat, method = "bonferroni")
pairwise <- pairwise.wilcox.test(density_df$density, density_df$MajorHabitat,
                     p.adjust.method = "BH")
pairwise

fig1 <- ggplot(density_df, aes(x = MajorHabitat, y = density)) +
  geom_boxplot(outlier.colour = "grey", outlier.shape = 16, outlier.size = 1) +
  geom_jitter(width = 0.3, size = 1, alpha = 0.2) +
  labs(y=bquote("Density (N per 100 m"^"2"*")")) +
  xlab("Habitat Type") +
  theme(axis.text.x=element_text(size=9))
fig1



#Figure 2: Density ~ Size class & Habitat
den_size <- select(density_df, -c(9:12))
den_size$density <- NULL
den_size <- den_size %>% 
  mutate_at(vars(2:8), list(~ . * 2))
library(reshape2)
melt_size <- melt(den_size, id = c("MajorHabitat"))
View(melt_size) 

size_sum <- ddply(melt_size, ~MajorHabitat*variable, summarise, 
                  mean = mean(value), 
                  sd = sd(value), 
                  n = length(value), 
                  SEM = sd(value)/sqrt(length(value)))

fig2.3 <- ggplot(size_sum, aes(x=variable, y= mean, fill = MajorHabitat))  +
  geom_col(position = position_dodge(),
           width = (0.75)) +
  geom_errorbar(aes(ymin = mean - SEM, ymax = mean + SEM), position = position_dodge(width = 0.9), width = 0.25) +
  labs(y=bquote("Average density (N per 100 m"^"2"*")"),
       fill = "Habitat Type") +
  xlab("Size Class (cm)") +
  scale_fill_manual(labels = c("Hardbottom", "Rocky reef", "Seagrass \n fringe"),
                    values=c("#F8766D", "#00B0F6", "#00BC59"))+
  scale_y_continuous(expand = c(0,0),
                     limits=c(0,6.5))

RedoFig2.2 <- ggplot(melt_size, aes(x = variable, y =value, fill = MajorHabitat)) +
  geom_boxplot(outlier.colour = "grey", outlier.shape = 16, outlier.size = 1) +
  labs(y=bquote("Density (N per 100 m"^"2"*")")) +
  xlab("Habitat type") +
  theme(axis.text.x=element_text(size=9)) +
  ylim(0, 13)
RedoFig2.2

fig2.3 <- ggplot(melt_size, aes(x=variable, y= mean, fill = MajorHabitat))  +
  geom_col(position = "fill",
           width = (0.75)) +
  labs(y=bquote("Proportion Density (N per 100 m"^"2"*")"),
       fill = "Habitat Type") +
  xlab("Size Class (cm)")
fig2.3

#Prelim Analysis
melt_prelim <- melt(prelim_overall, id = c("Habitat"))

prelim_fig <- ggplot(melt_prelim, aes(x=Habitat, y= value, fill = variable))  +
  geom_col(position = position_dodge(),
           width = (0.75)) +
  labs(y=bquote("Abundance"),
       fill = "Size Class (cm)") +
  xlab("Habitat Type")
prelim_fig
prelim_fig <- ggplot(melt_prelim, aes(x=Habitat, y= value, fill = variable))  +
  geom_col(position = position_dodge(),
           width = (0.75)) +
  labs(y=bquote("Abundance"),
       fill = "Size Class (cm)") +
  xlab("Habitat Type") +
  scale_fill_brewer(palette = "Dark2")
prelim_fig

redtail_df <- Overall_BeltTransectDF %>%
  filter(`Scientific name` == "sparisoma chrysopterum") %>%
  mutate(totalabun = rowSums(select(., 10:18))) %>%
  mutate(density = (totalabun*2)) 


kruskal.test(density ~ MajorHabitat, data = redtail_df)
dunn.test(density_df$density, g = density_df$MajorHabitat, method = "bonferroni")
pairwise.wilcox.test(density_df$density, density_df$MajorHabitat,
                                 p.adjust.method = "BH")







####################NBR: BENTHIC SUBSTRATES###############################
##########################################################################
modelpois <- glm(Abundance ~ MacaUNSP + Dictyota + TURF + SEAGRASS + CORAL + SPONGE, family = "poisson", data = glm_df)
summary(modelpois)
library(DHARMa)
resp <- simulateResiduals(modelpois, refit = T)
testDispersion(resp, plot = F)
plot(resp)

m1 <- glm.nb(Abundance ~ MacaUNSP + Dictyota + TURF + SEAGRASS + CORAL + SPONGE , data = glm_df)
summary(m1)
m2 <- update(m1, . ~ . - Dictyota)
anova(m1, m2)
m3 <- glm(Abundance ~ MacaUNSP + Dictyota + TURF + SEAGRASS + CORAL + SPONGE, family = "poisson", data = glm_df)
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)
qchisq(p=5.866e-69, df=8, lower.tail=FALSE)
#Chi-square 341.5115 which means NBR model is better

library(lmtest)
lrtest(m1, m3)
#I don't think this works for nbr models?
epiDisplay::poisgof(m1)
#I don't think this works for nbr models?


dat.resid <- sum(resid(m1, type = "pearson")^2)
dat.resid
m1$df.resid
1 - pchisq(dat.resid, m1$df.resid)

# Deviance (G2) residuals:
model$deviance
model$df.resid
1-pchisq(model$deviance, model$df.resid)

install.packages("pscl")
library(MASS)
negbin = glm.nb(Y~1)

# Reference website: https://www.uio.no/studier/emner/matnat/math/STK3100/h23/undervisningsmateriale/exercise_solutions/exercise_7.31.pdf
# Another: https://bookdown.org/ks6017/GLM_bookdown3/chapter-4-poisson-regression-and-extensions.html#poisson-distribution
library(pscl)
odTest(m1)
#Chi-square Test Statistic = 308.092, p < 0.0001
#Overdispersion is estimated in NBR 
install.packages("AER")
library(AER)
dispersiontest(m3)
#Dispersion 3.352

#INCLUDING SITE IN MODEL
install.packages("lme4")
trial <- glmer.nb(Density ~ MacaUNSP + Dictyota + TURF + SEAGRASS + CORAL + SPONGE +(1|SiteCode), family = binomial, data = glm_df)
summary(trial)
#Varience of random effect: 10.7


#NBR: Benthic Characteristics -
#####RUN ON DENSITY (N PER 100M2)
glm_df <- glm_df %>% mutate(Density = Abundance * 2)
trial <- glmer.nb(Density ~ MacaUNSP + Dictyota + TURF + SEAGRASS + CORAL + SPONGE +(1|SiteCode), family = binomial, data = glm_df)
m1 <- glm.nb(Density ~ MacaUNSP + Dictyota + TURF + SEAGRASS + CORAL + SPONGE , data = glm_df)
summary(m1)
m2 <- update(m1, . ~ . - Dictyota)
anova(m1, m2)
m3 <- glm(Density ~ MacaUNSP + Dictyota + TURF + SEAGRASS + CORAL + SPONGE, family = "poisson", data = glm_df)
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)
qchisq(p=7.971588e-239, df=8, lower.tail=FALSE)
#Chi square 1130
trial <- glmer.nb(Density ~ MacaUNSP + Dictyota + TURF + SEAGRASS + CORAL + SPONGE +(1|SiteCode), family = binomial, data = glm_df)
summary(trial)
plot_model(trial, type = "pred", terms=c("Dictyota"))
odTest(m1)
#Chi-square Test Statistic = 308.092, p < 0.0001
r.squaredGLMM(trial)

#Overdispersion is estimated in NBR 
#1089.0374 p-value = < 2.2e-16 





######GRAPH OF ABUNDANCE BY SITE###########
check <- abun_site %>%
  group_by(Site, MajorHabitat) %>%
  summarise(sum = mean(TotalAbundance))

plotsg <- abun_site %>%
  filter(MajorHabitat %in% c("Seagrass")) %>%
  ggplot(aes(x = Site, y = TotalAbundance)) +
  geom_bar(stat = "identity") +
  labs(x = "Site", y = "Density") +
  ggtitle("Seagrass") +
  theme_minimal()
plotsg

plothb <- abun_site %>%
  filter(MajorHabitat %in% c("Hardbottom")) %>%
  ggplot(aes(x = Site, y = TotalAbundance)) +
  geom_bar(stat = "identity") +
  labs(x = "Site", y = "Density") +
  ggtitle("Hardbottom") +
  theme_minimal()
plothb

plotrr <- abun_site %>%
  filter(MajorHabitat %in% c("Rocky Reef")) %>%
  ggplot(aes(x = Site, y = TotalAbundance)) +
  geom_bar(stat = "identity") +
  labs(x = "Site", y = "Density") +
  ggtitle("Rocky Reef") +
  theme_minimal()
plotrr

plotsg / plothb / plotrr 



####Benthic Cover Percent####

cpc_long <- pivot_longer(cpc_data,
                        cols = starts_with ("TRAN"),
                        names_to = "Category",
                        values_to = "Percent_Cover") %>%
                      filter(!is.na(Percent_Cover)) 

cpc_sum <- cpc_long %>%
  group_by(Habitat, COVER) %>%
  summarise(avg = mean(Percent_Cover))


plot_percentcov <- cpc_sum %>%
  ggplot(aes(x = Habitat, y = avg, fill = COVER)) +
  geom_col(position = "fill",
           width = (0.75)) +
  labs(x = "Site", y = "Density") +
  ggtitle("Seagrass") +
  theme_minimal()
plot_percentcov
#Need to figure out how we are combining this data#

####Chi Square####
chi <- density_df %>%
  group_by(Habitat) %>%l
  summarise(sum = mean(Percent_Cover))

contingency_table <- table(chi_df$`Size Class`, chi_df$Hardbottom)
contingency_table

column_sum <- sum(density_df$density)
#0-3
habitat <- c("Hardbottom", "Rocky Reef", "Seagrass")
size_class <- c("0-3")
count <- c(148, 140 , 108)
contingency_table <- matrix(count, nrow = length(habitat), ncol = length(size_class))
rownames(contingency_table) <- habitat
colnames(contingency_table) <- size_class
print(contingency_table)
size0_3 <- chisq.test(contingency_table)
size0_3
#p: 0.03

#4-5
habitat <- c("Hardbottom", "Rocky Reef", "Seagrass")
size_class <- c("4-5")
count <- c(550, 347 , 236)
contingency_table <- matrix(count, nrow = length(habitat), ncol = length(size_class))
rownames(contingency_table) <- habitat
colnames(contingency_table) <- size_class
size4_5 <- chisq.test(contingency_table)
size4_5
#p: < 0.0001

#6-10
habitat <- c("Hardbottom", "Rocky Reef", "Seagrass")
size_class <- c("4-5")
count <- c(606, 390, 250)
contingency_table <- matrix(count, nrow = length(habitat), ncol = length(size_class))
rownames(contingency_table) <- habitat
colnames(contingency_table) <- size_class
size6_10 <- chisq.test(contingency_table)
size6_10
#p: < 0.0001

#11-15
habitat <- c("Hardbottom", "Rocky Reef", "Seagrass")
size_class <- c("11-15")
count <- c(176, 172, 46)
contingency_table <- matrix(count, nrow = length(habitat), ncol = length(size_class))
rownames(contingency_table) <- habitat
colnames(contingency_table) <- size_class
chisq.test(contingency_table)
#p: < 0.0001

#16-20
habitat <- c("Hardbottom", "Rocky Reef", "Seagrass")
size_class <- c("16-20")
count <- c(56, 94, 34)
contingency_table <- matrix(count, nrow = length(habitat), ncol = length(size_class))
rownames(contingency_table) <- habitat
colnames(contingency_table) <- size_class
size16_20 <- chisq.test(contingency_table)
size16_20
#p: < 0.0001

#21-30
habitat <- c("Hardbottom", "Rocky Reef", "Seagrass")
size_class <- c("16-20")
count <- c(14, 24, 4)
contingency_table <- matrix(count, nrow = length(habitat), ncol = length(size_class))
rownames(contingency_table) <- habitat
colnames(contingency_table) <- size_class
size21_30<- chisq.test(contingency_table)
size21_30
#p: 0.0008

#31-40
habitat <- c("Hardbottom", "Rocky Reef", "Seagrass")
size_class <- c("16-20")
count <- c(0, 4, 4)
contingency_table <- matrix(count, nrow = length(habitat), ncol = length(size_class))
rownames(contingency_table) <- habitat
colnames(contingency_table) <- size_class
chisq.test(contingency_table)


long_data <- pivot_longer(chi_df, 
                          cols = c(Hardbottom, `Rocky Reef`, Seagrass),
                          names_to = "Category",
                          values_to = "Value")

chi_fig <- ggplot(long_data, aes(x= size, y= Value, fill = Category))  +
  geom_col(position = position_dodge(),
           width = (0.75)) +
  labs(y=bquote("Density (N per 100 m"^"2"*")"),
       fill = "Habitat Type") +
  xlab("Size Class (cm)")
chi_fig
long_data$size <- factor(long_data$size,levels = c("0_3", "4_5", "6_10", "11_15", "16_20", "21_30", "31_40"))




####NBR: DEPTH AND COMPLEXTY###
noout <- subset(nbr_complexdepth,Complexity < 100)
modelpois <- glm(Density ~ Complexity + AvgDepth, family = "poisson", data = noout)
summary(modelpois)
library(DHARMa)
resp <- simulateResiduals(modelpois, refit = T)
testDispersion(resp, plot = F)
plot(resp)
#SHOWS DISPERSION

m1 <- glm.nb(Density ~ Complexity + AvgDepth, data = noout)
summary(m1)
m2 <- update(m1, . ~ . - Complexity)
anova(m1, m2)
m1 <- glm.nb(Density ~ Complexity + AvgDepth, data = noout)
m3 <- glm(Density ~ Complexity + AvgDepth, family = "poisson", data = noout)
install.packages("AER")
library(AER)
dispersiontest(modelpois)
#Dispersion 3.441745 
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)
qchisq(p=5.16e-84, df=4, lower.tail=FALSE)
#Chi-square test 394.2 - ABOVE CRITICAL VALUE MEANS GOOD TEST
qchisq(0.05, df=4, lower.tail=FALSE)
summary(m1)

library(pscl)
odTest(m1)
#Critical value of test statistic at the alpha= 0.05 level: 2.7055 
#Chi-Square Test Statistic =  347.7721 p-value = < 2.2e-16 

m2 <- glm.nb(size1 ~ Complexity , data = noout)
colnames(lrdf)[1] = "size1"
s2 <- summary(m2)
s2
plot_model(m2, type = "pred", terms = "complexity", axis.title = "Density Prediction", title = "")
#0-3 not significant
size1 <- glmer.nb(size1 ~ Complexity + (1|SiteCode), family = binomial, data = noout)
summary(size1)
odTest(m2)
ggplot(chi_df, x = Habitat, y )

m3 <- glm.nb(size2 ~ Complexity , data = noout)
colnames(lrdf)[2] = "size2"
s3 <- summary(m3)
s3
plot_model(m3, type = "pred", terms = "complexity", axis.title = "Density Prediction", title = "")
#4-5 not significant
odTest(m3)
size2 <- glmer.nb(size2 ~ Complexity + (1|SiteCode), family = binomial, data = noout)
summary(size2)

m4 <- glm.nb(size3 ~ Complexity , data = noout)
colnames(lrdf)[3] = "size3"
s4 <- summary(m4)
s4
plot_model(m4, type = "pred", terms = "complexity", axis.title = "Density Prediction", title = "")
#6-10 significant
odTest(m4)
size3 <- glmer.nb(size3 ~ Complexity + (1|SiteCode), family = binomial, data = noout)
summary(size3)

m5 <- glm.nb(size4 ~ Complexity , data = noout)
s5 <- summary(m5)
s5
plot_model(m5, type = "pred", terms = "complexity", axis.title = "Density Prediction", title = "")
#11-15 significant
odTest(m5)
size4 <- glmer.nb(size4 ~ Complexity + (1|SiteCode), family = binomial, data = noout)
summary(size4)

m6 <- glm.nb(size5 ~ Complexity , data = noout)
s6 <- summary(m6)
s6
plot_model(m6, type = "pred", terms = "complexity", axis.title = "Density Prediction", title = "")
#16-20 significant
odTest(m6)
size5 <- glmer.nb(size5 ~ Complexity + (1|SiteCode), family = binomial, data = noout)
summary(size5)


m7 <- glm.nb(size6 ~ Complexity , data = noout)
s7 <- summary(m7)
s7
plot_model(m7, type = "pred", terms = "complexity", axis.title = "Density Prediction", title = "")
#21-30 significant but bad 
odTest(m7)
size6 <- glmer.nb(size6 ~ Complexity + (1|SiteCode), family = binomial, data = noout)
summary(size6)

# ANALYSIS NUMBER 6: DEPTH BY HABITAT TYPE
depth <- nbr_complexdepth 


sum_depth<- ddply(depth, ~MajorHabitat, summarise, 
                  mean = mean(DepthMeters),
                  sd = sd(DepthMeters), 
                  n = length(DepthMeters), 
                  SEM = sd(DepthMeters)/sqrt(length(DepthMeters)))
view(sum_depth)
kruskal.test(DepthMeters ~ MajorHabitat, data = depth)
dunnTest(DepthMeters ~ MajorHabitat, data = depth,
         method="bonferroni")
pairwise.wilcox.test(depth$DepthMeters, depth$MajorHabitat,
                     p.adjust.method = "BH")


depth_fig <- ggplot(depth, aes(x = MajorHabitat, y = DepthMeters)) +
  geom_boxplot(outlier.colour = "grey", outlier.shape = 16, outlier.size = 1) +
  geom_jitter(width = 0.3, size = 1, alpha = 0.2) +
  labs(y=bquote("Depth (m)")) +
  xlab("Habitat Type") +
  theme(axis.text.x=element_text(size=9))
depth_fig

#COMPLEXITY BETWEEN HABITAT TYPES
complex_sum  <- ddply(nbr_complexdepth, ~MajorHabitat, summarise,
                      mean = mean(Complexity), 
                      sd = sd(Complexity), 
                      n = length(Complexity), 
                      SEM = sd(Complexity)/sqrt
                      (length(Complexity)))

kruskal.test(Complexity ~ MajorHabitat, data = nbr_complexdepth)
#Kruskal-Wallis chi-squared = 129.35, df = 2, p-value <0.001 
#Significant Different Between Complexity in the Major Habitat Types
dunnTest(Complexity ~ MajorHabitat,
         data=nbr_complexdepth,
         method="bonferroni")

pairwise.wilcox.test(nbr_complexdepth$Complexity, nbr_complexdepth$MajorHabitat, 
                     p.adjust.method = "BH")

sigs <- c("a", "b", "c")

fig_complex <- ggplot(nbr_complexdepth, aes(x = MajorHabitat, y = Complexity)) +
  geom_boxplot(outlier.colour = "grey", outlier.shape = 16, outlier.size = 1) +
  geom_jitter(width = 0.3, size = 1, alpha = 0.2) +
  labs(y=bquote("Topographic Complexity (cm)")) +
  xlab("Habitat Type") +
  theme(axis.text.x=element_text(size=9))
fig_complex

plot_model(m_depthcomplex, type = "pred", terms = "Complexity", axis.title = "Density Prediction", title = "")

#PROPORTION GRAPH OF BENTHIC SUBSTRATES

#cATEGORIES FOR GRAPH:
#SEAGRASS, OTHER, DICYTOTA, TURF ALGAE, CORAL, SPONGE, NON-LIVING, MACROALGAE UNSPECIFIED

#rowsum(GORGONIAN, ZOANTHID, CCA, CYANOBACTERIA, OTHER)
benthicdf$allother <- rowSums( benthicdf[, c(7, 9, 11, 12, 19) ] )

benthic$Habitat[benthic$Habitat == "Seagrass"] <- "Seagrass Fringe"

benthic <- benthicdf %>%
  subset(select = c(Habitat, CORAL, SPONGE, TURF, Nonliving, SEAGRASS, allother, MacaUNSP, DICT)) %>%
  rename("Coral" = "CORAL",
         "Sponge" = "SPONGE", 
         "Turf Algae" = "TURF", 
         "Non-living" = "Nonliving",
         "Other" = "allother",
         "Macroalgae Unsp" = "MacaUNSP", 
         "Seagrass" = "SEAGRASS",
         "Dictyota" = "DICT") %>%
  melt(id = ("Habitat")) %>%
  na.omit()
benthic$Habitat[benthic$Habitat == "Seagrass"] <- "Seagrass Fringe"


benthic_graph <-  ggplot(benthic, aes(x=Habitat, y= value, fill = variable))  +
  geom_col(position = "fill",
           width = (0.75)) +
  labs("Percent Cover", fill = "Biological Substrate") +
  xlab("Habitat Type") +
  ylab("Percent Cover") +
  theme(axis.title.y = element_text(size = 10))
benthic_graph

level_order_var <- c("Seagrass", "Other", "Sponge", "Coral", "Macroalgae Unsp", "Dictyota","Turf Algae","Non-living")
benthic$variable <- factor(benthic$variable, levels=level_order_var)


ggplot() +
  geom_point(data = noout, aes(x = Complexity,
                               y = Density,
  ),
  alpha = 0.2) +
  labs(y=bquote("Predicted density (N per 50 m"^"2"*")")) +
  labs(x = "Complexity (cm)") +
  geom_path(data = CI_comp, aes(x = x, y = predicted, color = "red", )) +
  geom_ribbon(data = CI_comp, aes(x = x, ymin = conf.low,
                                  ymax = conf.high),
              alpha = 0.5) 


ggplot() +
  geom_point(data = glm_df, aes(x = Dictyota,
                                y = Abundance,
  ),
  alpha = 0.2) +
  labs(y=bquote("Predicted Density (N per 50 m"^"2"*")")) +
  labs(x = "Percent Cover Dictyota Algae") +
  geom_path(data = CI_ben, aes(x = x, y = predicted, color = "red")) +
  geom_ribbon(data = CI_ben, aes(x = x, ymin = conf.low,
                                 ymax = conf.high),
              alpha = 0.5)
