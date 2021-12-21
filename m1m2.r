source('figures_source.r')
dat <- load_main('ldat')

# m1m2 --------------------------------------------------------------------

h1 <- read.csv('data_sc/escape/m1m2.csv') %>% select(-X)
h2 <- na.omit(h1)

# overall cell distributions
ggplot(h2,aes(clinical,m1m2))+geom_violin(aes(fill=clinical))#+facet_wrap(~cluster_fine)

# overall distribution of scores
ggplot(h2,aes(m1m2))+geom_density(adjust=0.01)
  # hard to define M1 and M2 subsets but int/cd16 are generally more inflammatory so checks out

# by patient
h3 <- group_by(h2,patient,timepoint,cluster_fine) %>% 
  summarize(m1m2=mean(m1m2)) %>% add_clinical()
group_by(h3,clinical,cluster_fine) %>% summarize(se=sd(m1m2)/sqrt(n()),m1m2=mean(m1m2)) %>%
  ggplot(aes(cluster_fine,m1m2,fill=clinical))+
  geom_col(position=position_dodge(width=1))+
  geom_errorbar(aes(ymin=m1m2-se,ymax=m1m2+se),width=0.2,lwd=1.2,
                position=position_dodge(width=1))+
  geom_point(data=h3,aes(cluster_fine,m1m2,color=clinical),
             position=position_dodge(width=1))+
  scale_color_manual(values=c('black','black','black','black'))

# timepoint split by patient
group_by(h2,patient,timepoint,cluster_fine) %>% 
  summarize(m1m2=mean(m1m2)) %>% add_clinical() %>% filter(clinical!='Healthy') %>% 
  ggplot(aes(timepoint,m1m2,fill=clinical))+
  stat_summary(fun=mean, geom = "bar",position=position_dodge(width=1)) +
  stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position=position_dodge(width=1))+
  facet_wrap(~cluster_fine)

# timepoint propotional change split by patient
group_by(h2,patient,timepoint,cluster_fine) %>% 
  summarize(m1m2=mean(m1m2)) %>% add_clinical() %>% filter(clinical!='Healthy') %>% 
  pivot_wider(id_cols=c(patient,cluster_fine),names_from=timepoint,values_from=m1m2) %>% 
  mutate(d13=(`3`-`1`)/`1`,d15=(`5`-`1`)/`1`) %>% add_clinical() %>% 
  # ggplot(aes(cluster_fine,d13,fill=clinical))+
  ggplot(aes(cluster_fine,d15,fill=clinical))+
  geom_hline(yintercept=0)+geom_boxplot()
  

# NK cell counts ----------------------------------------------------------------------

h1 <- (dat$timepoint %in% c('1','3','5'))&
  (interaction(dat$patient,dat$timepoint) %!in% c('20.5','23.1'))&
  ((dat$cluster_coarse %!in% c('Other','Junk'))|(dat$cluster_fine=='Hematopoetic Progenitors'))
h1 <- table(dat$cluster_fine[h1],
      interaction(dat$timepoint[h1],dat$patient[h1],drop=T)) %>% 
  as.data.frame() %>% separate(Var2,c('timepoint','patient')) %>% 
  rename(type=Var1) %>% 
  group_by(patient,timepoint) %>% mutate(Freq=proportions(Freq)) %>% 
  mutate(timepoint=factor(if_else(timepoint=='1','pre','post'),
                          levels=c('pre','post'))) %>% 
  pivot_wider(id_cols=c(type,patient),names_from=timepoint,values_from=Freq,
              values_fn=mean) %>% 
  pivot_longer(c(pre,post),names_to='timepoint',values_to='freq') %>% 
  add_clinical() %>% filter(type=='CD16 NK')
ggplot(h1,aes(factor(timepoint,c('pre','post')),freq,color=clinical,
              group=patient))+geom_path()+xlab('')+
  theme_prism()+scale_color_manual(values=cols_clinical)

pivot_wider(h1,patient,timepoint,values_from=freq) %>% na.omit() %>% 
  mutate(delta=(post-pre)/pre) %>% mutate(timepoint='post') %>% 
  right_join(h1,by=c('timepoint','patient')) %>% 
  mutate(delta=if_else(timepoint=='pre',0,delta)) %>% 
  ggplot(aes(factor(timepoint,c('pre','post')),delta,color=clinical,
                group=patient))+geom_path()+
  theme_prism()+scale_color_manual(values=cols_clinical)+facet_wrap(~clinical)+
  coord_cartesian(ylim=c(-0.5,0.5))


# macro stuff -------------------------------------------------------------

# escape GSEA enrich scores for these sets for all monocytes in ldat
mac_diff <- readRDS('data_sc/escape/GOBP_MACROPHAGE_DIFFERENTIATION.rds') %>% 
  rename('mac_diff'=val)
mono_chemo <- readRDS('data_sc/escape/GOBP_POSITIVE_REGULATION_OF_MONOCYTE_CHEMOTAXIS.rds') %>% 
  rename('mono_chemo'=val)
mono_diff <- readRDS('data_sc/escape/GOBP_POSITIVE_REGULATION_OF_MONOCYTE_DIFFERENTIATION.rds') %>% 
  rename('mono_diff'=val)

h1 <- full_join(mac_diff,mono_chemo,by='barcode') %>% 
  full_join(mono_diff,by='barcode') %>% na.omit() %>% 
  left_join(select(dat,barcode,cluster_fine,patient,timepoint),by='barcode') %>% 
  na.omit() %>% select(-barcode) %>% add_clinical() %>% 
  pivot_longer(c(mac_diff,mono_chemo,mono_diff))

# overall cell distributions
ggplot(h1,aes(clinical,value)) %>% prism()+geom_violin(aes(fill=clinical))+
  facet_wrap(~name)+theme(axis.text.x=element_blank())

# overall distribution of scores
ggplot(h1,aes(value,color=name))+geom_density(adjust=0.5)

#### by patient
group_by(filter(h1,name=='mac_diff'),patient,timepoint,cluster_fine) %>% 
  summarize(value=mean(value)) %>% add_clinical() %>% 
  group_by(clinical,cluster_fine) %>% summarize(se=sd(value)/sqrt(n()),value=mean(value)) %>% 
  ggplot(aes(cluster_fine,value,fill=clinical)) %>% prism()+
  geom_col(position=position_dodge(width=1))+
  geom_errorbar(aes(ymin=value-se,ymax=value+se),width=0.2,lwd=1.2,
                position=position_dodge(width=1))

# timepoint split by patient
group_by(filter(h1,name=='mac_diff'),patient,timepoint,cluster_fine) %>% 
  summarize(value=mean(value)) %>% add_clinical() %>% filter(clinical!='Healthy') %>% 
  ggplot(aes(timepoint,value,fill=clinical)) %>% prism()+
  stat_summary(fun=mean, geom = "bar",position=position_dodge(width=1)) +
  stat_summary(fun.data = mean_se, geom = "errorbar",width=0.2,position=position_dodge(width=1))+
  facet_wrap(~cluster_fine)

# timepoint propotional change split by patient
group_by(filter(h1,name=='mac_diff'),patient,timepoint,cluster_fine) %>% 
  summarize(value=mean(value)) %>% add_clinical() %>% filter(clinical!='Healthy') %>% 
  pivot_wider(id_cols=c(patient,cluster_fine),names_from=timepoint,values_from=value) %>% 
  mutate(d13=(`3`-`1`)/`1`,d15=(`5`-`1`)/`1`) %>% add_clinical() %>% 
  # ggplot(aes(cluster_fine,d13,fill=clinical))+
  ggplot(aes(cluster_fine,d15,fill=clinical)) %>% prism()+
  geom_hline(yintercept=0)+geom_boxplot()+
  ylab('Prop change timepoint 1-5')


# NK stuff ----------------------------------------------------------------

# escape GSEA enrich scores for these sets for all monocytes in ldat
h1 <- cbind(
  readRDS('data_sc/escape/GOBP_POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_ACTIVATION.rds') %>% rename(activation=val),
  readRDS('data_sc/escape/GOBP_POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_CHEMOTAXIS.rds') %>% select(chemotaxis=val),
  readRDS('data_sc/escape/GOBP_POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY.rds') %>% select(cytotoxic=val),
  readRDS('data_sc/escape/GOBP_POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL.rds') %>% select(tumor_response=val),
  readRDS('data_sc/escape/GOBP_POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_PROLIFERATION.rds') %>% select(proliferation=val)
) %>% left_join(select(dat,patient,timepoint,clinical,cluster_fine,barcode),by='barcode') %>% 
  na.omit() %>% select(-barcode) %>% 
  pivot_longer(-c(clinical,patient,timepoint,cluster_fine),names_to='pathway')

#### see m1m2_nk.rmd for figures and viz