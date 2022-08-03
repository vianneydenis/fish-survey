library(openxlsx)
library(dplyr)
library(rfishbase)
library(reshape)
library(ggplot2)
library(ggpattern)
library(ggrepel)
library(ggpubr)
library(geodist)
library(vegan)
library(betapart)
library(VennDiagram)
library(NbClust)
library(corrplot)
library(factoextra)
library(iNEXT)
library(stringr)
library(fishtree)
library(ggtree)
library(caper)
library(picante)
library(pairwiseAdonis)
library(rstatix)
library(nlme)
library(multcomp)
library(PMCMRplus)



rm(list=ls())
NE = read.csv("Data/Site.csv") # site information
sp_D = read.csv("Data/DOV.csv") # diver-operated videos (DOV) data
sp_U = read.csv("Data/UVC.csv") # underwater visual census (UVC) data
sp_eDNA = read.xlsx("Data/eDNA.xlsx") # Environmental DNA (eDNA) data
bc_mt = read.csv("Data/BV.csv", row.names = 1) # Benthic variables
sp_lev = read.csv("Data/level.csv") # Species' vertical position in the water column 
sp_lwab = read.csv("Data/Lw_ab.csv") # Biomass coefficients
flow = read.csv("Data/flow.csv") # Water flow data
gd = read.csv("Data/Pairwise geographic distance.csv", header = T, row.names = 1, check.names = F) ## Pairwise geographic distance (avoiding the land)

sp_sv_D = validate_names(sp_D$Species)
sp_sv_U = validate_names(sp_U$Species)
sp_sv_e = validate_names(sp_eDNA[, 5])

abun_D = as.data.frame(cast(sp_D[,-1], Site ~ Species,
                          value='Number', fun.aggregate = sum))
rownames(abun_D) = NE$Code
abun_D = abun_D[,-1]
rich_D = abun_D
rich_D[rich_D>0] = 1

abun_D = abun_D[,-1]
abun_U = as.data.frame(cast(sp_U[,-1], Site ~ Species,
                            value='Number', fun.aggregate = sum))
rownames(abun_U) = NE$Code
abun_U = abun_U[,-1]
rich_U = abun_U
rich_U[rich_U>0] = 1

sp_eDNA = sp_eDNA[sp_eDNA$Ratio>0.0001,] # remove reads with coverage < 0.01%
abun_e = as.data.frame(cast(sp_eDNA, Site ~ Valid_as, value='Ratio', fun.aggregate = sum))
row.names(abun_e) = abun_e$Site
abun_e = abun_e[,-1]
rich_e = abun_e
rich_e[rich_e>0] = 1

rich_3 = full_join(full_join(rich_D, rich_U), rich_e)
rich_3 = cbind(data.frame(Method = rep(c("DOV", "UVC", "eDNA"), each = 21),
                    Site = rep(NE$Code, 3)), rich_3)
rich_3[is.na(rich_3)] = 0

sp_bio = rbind(sp_D, sp_U)
# retrieve coefficient a
sp_bio$a = NA
for (n in 1:nrow(sp_bio)) {
  sp_bio$a[n] = sp_lwab$mean_a[which(sp_bio$Species[n] == sp_lwab$Species)]
}
# retrieve coefficient b
sp_bio$b = NA
for (n in 1:nrow(sp_bio)) {
  sp_bio$b[n] = sp_lwab$mean_b[which(sp_bio$Species[n] == sp_lwab$Species)]
}
# calculate biomass (g/m^2) based on total length-weight equation (g = a*TL^b)
sp_bio$biomass = (sp_bio$a*(sp_bio$Length)^sp_bio$b)/250
sp_bio = na.omit(sp_bio)

ncol(rich_D) # DOV: 83 spp
ncol(rich_U) # UVC: 111 spp
ncol(rich_e) # eDNA: 383 spp

m = c(1:42) # set the min and max sampling effort
# generate species interpolation and extrapolation curves based on occurrence data
out = iNEXT(list(DOV = t(rich_D), 
                 UVC = t(rich_U), 
                 eDNA = t(rich_e)), 
            q = 0, datatype = "incidence_raw", size = m)
names(out$iNextEst) = c("DOV", "UVC", "eDNA") # label each method

# find out the number of site when reaching the saturation (<1% increase)
d0 = out$iNextEst$DOV$qD[-42] # 1-41
d1 = out$iNextEst$DOV$qD[-1] # 2-42
sp_d = min(which(((d1-d0)/d0*100)<1))+1 # 39 sites

u0 = out$iNextEst$UVC$qD[-42] # 1-41
u1 = out$iNextEst$UVC$qD[-1] # 2-42
sp_u = min(which(((u1-u0)/u0*100)<1))+1 # 40 sites

e0 = out$iNextEst$eDNA$qD[-42]
e1 = out$iNextEst$eDNA$qD[-1]
sp_e = min(which(((e1-e0)/e0*100)<1))+1 # 34 sites

dp = data.frame(
  Value = c(ncol(rich_D), ncol(rich_U), ncol(rich_e),
            round(out$iNextEst$DOV$qD[sp_d]),
            round(out$iNextEst$UVC$qD[sp_u]),
            round(out$iNextEst$eDNA$qD[sp_e])),
  X = c(rep(21,3), sp_d, sp_u, sp_e),
  Y = c(ncol(rich_D)-40, ncol(rich_U)+40, ncol(rich_e)+40,
        round(out$iNextEst$DOV$qD[sp_d]-40),
        round(out$iNextEst$UVC$qD[sp_u]+40),
        round(out$iNextEst$eDNA$qD[sp_e])+60),
  Method = as.factor(rep(c("DOV", "UVC", "eDNA"), 2)),
  Type = rep(c("Observed", "Saturated"), each = 3)
)
dp$Text = str_c("(", dp$X, ", ", dp$Value, ")", sep = "")

gg = ggiNEXT(out, type=1)+
  theme(plot.title = element_text(face = "bold", size = (15), hjust = 0.5),
        axis.text = element_text(colour = "black", face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.position = "right", legend.key=element_blank())+
  labs(x = "Sampling units", y = "Species richness")
gg$layers[[1]] = NULL

# Figure 2. Species rarefaction curves.
gg+ geom_point(data = dp, aes(x = X, y= Value, col = Method, shape = Type), size = 5)+
  geom_text(data = dp, aes(x = X, y= Y, label = Text, col = Method), size = 5)+
  scale_shape_manual(values = 16:17)+
  scale_x_continuous(breaks = c(0, 21, 42))

# extract species list for each method
venn_e = colnames(rich_e) # eDNA: 383 spp
venn_D = colnames(rich_D) # DOC: 83 spp
venn_U = colnames(rich_U) # UVC: 111 spp

venn_object = venn.diagram(list(venn_e, venn_U, venn_D), filename = NULL,
                            category.names = c("eDNA", "UVC", "DOV"),
                            output=TRUE,
                            # Output features
                            imagetype="svg" ,
                            height = 100 , 
                            width = 100 , 
                            resolution = 300,
                            compression = "lzw",
                            # Circles
                            lwd = 2,
                            lty = 'blank',
                            fill = c("#7CAE00", "#00BFC4", "#F8766D"),
                            # Numbers
                            cex = 1,
                            fontface = "bold",
                            fontfamily = "sans",
                            cat.cex = 1,
                            cat.col = c("#7CAE00", "#00BFC4", "#F8766D"),
                            cat.fontface = "bold",
                            cat.pos = c(-30, 30, 180),
                            cat.dist = c(0.08, 0.08, 0.055),
                            cat.fontfamily = "sans")
# Figure 3. Number of species detected by each method
grid.newpage()
grid.draw(venn_object)

# include information in richness matrix
lev_D = table(left_join(data.frame(Species = colnames(rich_D)), sp_lev, "Species")$level)
lev_U = table(left_join(data.frame(Species = colnames(rich_U)), sp_lev, "Species")$level)
lev_e = table(left_join(data.frame(Species = colnames(rich_e)), sp_lev, "Species")$level)
lev_3 = data.frame(DOV = c(round(lev_D/sum(lev_D)*100), 0), 
                   UVC = c(round(lev_U/sum(lev_U)*100), 0), 
                   eDNA = as.vector(round(lev_e/sum(lev_e)*100)),
                   row.names = c("Benthic", "Bentho-pelagic", "Pelagic"))
# prepare information for plotting
df = melt(lev_3, variable_name = "Method")
df$Level = rep(c("Benthic", "Bentho-pelagic", "Pelagic"), 3)
# Figure S1. Composition of fishes’ level in the water column among methods.
ggplot(df[-c(3,6),], aes(y = value, x = Method, fill = Level)) +
  geom_bar(stat="identity", width= 0.6)+
  labs(y = "Percentage (%)", x = "Method")+
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        axis.title = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.position = "bottom", legend.key=element_blank())+
  geom_text(aes(label = value), size = 4, position =  position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("#fdae61", "#a6cee3", "#1f78b4"))

# retrieve phylogentic tree from Fish Tree of Life
tre = fishtree_phylogeny(colnames(rich_3[-1:-2]), type = "phylogram") # only 335 species were available

# Screen out species detected by each method
n3 = list(colnames(rich_D), colnames(rich_U), colnames(rich_e))
ss = vector(mode = "list", length = 3)
for (i in 1:3) {
  ss[[i]] = strsplit(n3[[i]], split = " ") # split genus and species
  # change species name into Genus_species
  for (j in 1:length(n3[[i]])) {
    n3[[i]][j] = stringr::str_c(ss[[i]][[j]][1], ss[[i]][[j]][2], sep = "_")
  }
  n3[[i]] = intersect(n3[[i]], tre$tip.label)
  # DOV: 72 of 83 species were available
  # UVC: 100 of 111 species were available
  # eDNA: 294 of 413 species were available
  }
# Add grouping information to the tree
tre_DOV = groupOTU(tre, n3[[1]], group_name = "DOV")
tre_UVC = groupOTU(tre, n3[[2]], group_name = "UVC")
tre_eDNA = groupOTU(tre, n3[[3]], group_name = "eDNA") 

# Plot trees
gt_D = ggtree(tre_DOV, branch.length='none', layout='circular',aes(color = DOV))+
  geom_tiplab(size=1, aes(angle=angle))+
  scale_color_manual(values=c("black", "deeppink"))+
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", size = (24), color = "#F8766D", hjust = 0.5))+
  ggtitle("DOV")
gt_U = ggtree(tre_UVC, branch.length='none', layout='circular',aes(color = UVC))+
  geom_tiplab(size=1, aes(angle=angle))+
  scale_color_manual(values=c("black", "blue"))+
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", size = (24), color = "#619CFF", hjust = 0.5))+
  ggtitle("UVC")
gt_e = ggtree(tre_eDNA, branch.length='none', layout='circular',aes(color = eDNA))+
  geom_tiplab(size=1, aes(angle=angle))+
  scale_color_manual(values=c("black", "#7CAE00"))+
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", size = (24), color = "#00BA38", hjust = 0.5))+
  ggtitle("eDNA")

# Figure S2. Phylogenetic trees.
ggarrange(gt_D, gt_U, gt_e, ncol = 3, nrow = 1, common.legend = F,
          labels = c("(a)", "(b)", "(c)"))

# generate occurrence data (methods by species) 
PD_sp = aggregate(.~Method, rich_3[,-2], FUN = sum)
row.names(PD_sp) = PD_sp$Method
PD_sp = PD_sp[,-1]
PD_sp[PD_sp>0] = 1
ss = strsplit(colnames(PD_sp), split = " ") # split genus and species
# change species name into Genus_species
for (i in 1:ncol(PD_sp)) {
  colnames(PD_sp)[i] = stringr::str_c(ss[[i]][1], ss[[i]][2], sep = "_") 
}
# Standardized effect size of PD
SPD = ses.pd(PD_sp, tre, include.root = F, null.model = "taxa.labels")
SPD
# only PD in eDNA requires size standardization (p>0.05)


D_sp = as.data.frame(t(PD_sp))
D_sp$Species = row.names(D_sp)
Fish = comparative.data(tre, D_sp, Species) # combine phylogeny with occurence data


DOVPhyloD = phylo.d(Fish, binvar = DOV, permut = 999)
DOVPhyloD


UVCPhyloD = phylo.d(Fish, binvar = UVC, permut = 999) 
UVCPhyloD


eDNAPhyloD = phylo.d(Fish, binvar = eDNA, permut = 999) 
eDNAPhyloD
# Degree of disperse: eDNA > DOV > UVC

## nMDS on occrrence data across methods
nmds_3 = metaMDS(rich_3[, -1:-2], distance='jaccard',trymax=999)
nmds_3$stress # stress value = 0.1482547
# extract data from nMDS for plotting
data.scores_3 = as.data.frame(scores(nmds_3)$sites) # extract nMDS value of the 21 sites
data.scores_3$Method = rich_3$Method
data.scores_3$Site = rich_3$Site
data.scores_3$Area = NE$Area
species.scores_3 = as.data.frame(scores(nmds_3,"species")) # extract nMDS value of all species
species.scores_3$species = row.names(species.scores_3)
# Figure 4. nMDS of fish assemblages among sites across methods.
ggplot(data.scores_3, aes(x = NMDS1, y = NMDS2))+
  theme(plot.title = element_blank(),
        axis.text = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right",
        axis.title = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())+
  geom_point(size = 4, aes(colour = Method, shape = Area))+
  geom_text_repel(aes(label = Site))+
  labs(x = "NMDS1", y = "NMDS2", colour = "Method")+
  annotate(geom="text", x=0.75, y=0.85, label= paste("Stress =", round(nmds_3$stress,2)),
           size = 6, fontface = "bold")+
  scale_shape_discrete(labels = c("Northern coast", "Outlying islands"))

# PERMANOVA
PERM = adonis(rich_3[,-1:-2]~ rich_3$Method, strata = rich_3$Site, 
              permutations = 999, method="jaccard")
PERM$aov.tab # significant difference among methods
# Dispersion test
rich_3_j = vegdist(rich_3[, -1:-2], method="jaccard") 
mod = betadisper(rich_3_j, rich_3$Method, type=c('centroid'))
permutest (mod, pairwise=T) # p>0.05: no dispersion among methods
# Centroid test
pairwise.adonis(rich_3[, -1:-2], factors = rich_3$Method, sim.method = "jaccard") # eDNA-DOV, p<0.001; eDNA-UVC, p<0.001; DOV-UVC, p>0.05

# PERMANOVA
PERM_D = adonis(rich_D~ NE$Area, permutations = 999, method="jaccard")
PERM_D$aov.tab # DOV: significant difference between areas; p<0.001
# Dispersion test
rich_D_j = vegdist(rich_D, method="jaccard") 
mod = betadisper(rich_D_j, NE$Area, type=c('centroid'))
permutest (mod, pairwise=T) # p>0.05: no dispersion between areas
# Centroid test
pairwise.adonis(rich_D, factors = NE$Area, sim.method = "jaccard") # p<0.001; significant difference in centroids between areas

# PERMANOVA
PERM_U = adonis(rich_U~ NE$Area, permutations = 999, method="jaccard")
PERM_U$aov.tab # UVC: significant difference between areas; p<0.001
# Dispersion test
rich_U_j = vegdist(rich_U, method="jaccard") 
mod = betadisper(rich_U_j, NE$Area, type=c('centroid'))
permutest (mod, pairwise=T) # p>0.05: no dispersion between areas
# Centroid test
pairwise.adonis(rich_U, factors = NE$Area, sim.method = "jaccard") # p<0.001; significant difference in centroids between areas


PERM_e = adonis(rich_e~ NE$Area, permutations = 999, method="jaccard")
PERM_e$aov.tab # eDNA: significant difference between areas; p<0.01
rich_e_j = vegdist(rich_e, method="jaccard") 
mod = betadisper(rich_e_j, NE$Area, type=c('centroid'))
permutest (mod, pairwise=T) #p<0.001: there is dispersion between areas
pairwise.adonis(rich_e, factors = NE$Area, sim.method = "jaccard") # p<0.01; significant difference in centroids between areas

beta_D = beta.multi(rich_D, index.family="jaccard") # DOV
beta_D$beta.JAC # overall, 0.9433685
beta_D$beta.JTU # turnover, 0.9212947
beta_D$beta.JNE # nestedness, 0.02207375

beta_U = beta.multi(rich_U, index.family="jaccard") # UVC
beta_U$beta.JAC # overall, 0.944301
beta_U$beta.JTU # turnover, 0.929008
beta_U$beta.JNE # nestedness, 0.015293


beta_e = beta.multi(rich_e, index.family="jaccard") # eDNA
beta_e$beta.JAC # overall, 0.9459638
beta_e$beta.JTU # turnover, 0.9285283
beta_e$beta.JNE # nestedness, 0.01743549


# preparing plotting information
beta_plot = data.frame(Component = rep(c("turnover", "nestedness"), 3), 
                       Method = rep(c("DOV", "UVC", "eDNA"), each = 2), 
                       Beta_diversity = round(c(beta_D$beta.JTU, beta_D$beta.JNE, beta_U$beta.JTU, beta_U$beta.JNE, beta_e$beta.JTU, beta_e$beta.JNE), 2))
beta_plot$Component = factor(beta_plot$Component, levels=unique(beta_plot$Component))

s = ggplot(beta_plot)+ 
  aes(x = Method, y = Beta_diversity, fill = Method, pattern = Component)+
  theme(plot.title = element_text(face = "bold", size = (15), hjust = 0.5),
        axis.title = element_text(face = "bold", size = 14, colour = "black"),
        axis.text = element_text(colour = "black", face = "bold", size = 12),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),         legend.position = "right",
        legend.key=element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
  coord_cartesian(ylim=c(0,1))+
  geom_bar_pattern(stat="identity", width= 0.6, color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6)+
  scale_pattern_manual(values = c(nestedness = "stripe", turnover = "none"))+
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))+
  geom_text(aes(label = Beta_diversity), size = 3, hjust = 0.5, vjust = -0.5, position = "stack")+
  labs(x = "Method", y = "Beta diversity")
# Figure S3. The spatial patterns of β-diversity among methods.
s+guides(fill = "none")

beta_p_D = beta.pair(rich_D, index.family="jaccard") # DOV
beta_p_U = beta.pair(rich_U, index.family="jaccard") # UVC
beta_p_e = beta.pair(rich_e, index.family="jaccard") # eDNA

# pairwise beta diversity ~ pairwise geographic distance
mod = list(mantel(gd,as.matrix(beta_p_D[[3]])),
           mantel(gd,as.matrix(beta_p_D[[1]])),
           mantel(gd,as.matrix(beta_p_D[[2]])),
           mantel(gd,as.matrix(beta_p_U[[3]])),
           mantel(gd,as.matrix(beta_p_U[[1]])),
           mantel(gd,as.matrix(beta_p_U[[2]])),
           mantel(gd,as.matrix(beta_p_e[[3]])),
           mantel(gd,as.matrix(beta_p_e[[1]])),
           mantel(gd,as.matrix(beta_p_e[[2]])))

# generate empty lists
bd_D = vector(mode = "list", length = 3)
bd_U = vector(mode = "list", length = 3)
bd_e = vector(mode = "list", length = 3)
gd_vec = vector()
# collapsing pairwise data.frame into vector for plotting
for (n in 1:20) {
  gd_vec = c(gd_vec, gd[n, (n+1):21])
  for (m in 1:3) {
    bd_D[[m]] = c(bd_D[[m]], as.matrix(beta_p_D[[m]])[n, (n+1):21])
    bd_U[[m]] = c(bd_U[[m]], as.matrix(beta_p_U[[m]])[n, (n+1):21])
    bd_e[[m]] = c(bd_e[[m]], as.matrix(beta_p_e[[m]])[n, (n+1):21])
    }
}
gd_vec = as.numeric(gd_vec)

ds = list(bd_D[[3]], bd_D[[1]], bd_D[[2]],
          bd_U[[3]], bd_U[[1]], bd_U[[2]],
          bd_e[[3]], bd_e[[1]], bd_e[[2]])

z_lab = rep(c(0.1, 0.1, 0.9), 3) # position of mantel results on the plot
col3 = rep(c("#F8766D", "#619CFF", "#00BA38"), each = 3)
YLAB = rep(c("Pairwise β-diversity", "Pairwise turnover","Pairwise nestedness"), 3)
gg = list()
mt = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)")
for (n in 1:9) {
  df = data.frame(x = gd_vec, y = ds[[n]])
  gg[[n]] = ggplot(df, aes(x = x, y = y))+
    labs(x = "Geographic distance (km)", y = YLAB[n])+
    theme(plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = -0.15), 
          axis.title = element_text(face = "bold", size = 14, colour = "black"),
          axis.text = element_text(colour = "black", face = "bold", size = 12), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "bottom",
          legend.key=element_blank(), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))+
    ggtitle(mt[n])+
    scale_y_continuous(limits = c(0, 1))+
    geom_point(col = col3[n])+
    geom_smooth(formula = y~x, method = "lm")+
    geom_label(x = 65, y = z_lab[n], 
               label = paste("r =", round(mod[[n]][["statistic"]],2), ", p =" , round(mod[[n]][["signif"]],3)))
}
# Figure 5. Mantel correlations between pairwise β-diversity and the pairwise geographic distance.
ggarrange(gg[[1]], gg[[2]], gg[[3]], gg[[4]], gg[[5]],
          gg[[6]], gg[[7]], gg[[8]], gg[[9]],
          nrow = 3, ncol = 3, common.legend = F)

df = data.frame(Richness = rowSums(rich_3[,-1:-2]), Method = as.factor(rich_3$Method), Site = rep(NE$Code, 3))
# mean and sd
df %>% group_by(Method) %>% summarise(Mean = round(mean(Richness),2), SD = round(sd(Richness),2))
# One-way repeated measures ANOVA
mod = lme(Richness ~ Method, random = ~1|Site/Method, data=df)
anova(mod) # p<0.001
summary(glht(mod,linfct = mcp(Method = "Tukey"))) # p<0.001: eDNA-DOV, eDNA-UVC

# normalize benthic composition in each site as constrained variables
benthos_mt = bc_mt[,-11] %>% mutate_all(scale) %>% as.data.frame()
# calculate the variance inflation factors for each variable
diag(solve(cor(benthos_mt))) # VIF>10: CCA, hard coral, turf
# remove the variable, turf, with the largest VIF
diag(solve(cor(benthos_mt[,-9]))) # no multicollinearity
benthos_mt = benthos_mt[,-9]
# aggregate biomass data for each site and method
bioV = cast(Method+Site~Species, data = sp_bio, value='biomass', fun.aggregate = sum)

# occurrence
rich_D.jac = as.matrix(vegdist(rich_D, method = "jaccard"))
dbRDA.mat_DJ = capscale(rich_D.jac ~ ., benthos_mt, comm = rich_D)
# abundance
abun_D.b = as.matrix(vegdist(abun_D, method = "bray"))
dbRDA.mat_DB = capscale(abun_D.b ~ ., benthos_mt, comm = abun_D)
# biomass
bio_D.b = as.matrix(vegdist(bioV[1:21,-1:-2], method = "bray"))
dbRDA.mat_Db = capscale(bio_D.b ~ ., benthos_mt, comm = bioV[1:21,-1:-2])

# occurrence
rich_U.jac = as.matrix(vegdist(rich_U, method = "jaccard"))
dbRDA.mat_UJ = capscale(rich_U.jac ~ ., benthos_mt, comm = rich_U)
# abundance
abun_U.b = as.matrix(vegdist(abun_U, method = "bray"))
dbRDA.mat_UB = capscale(abun_U.b ~ ., benthos_mt, comm = abun_U)
# biomass
bio_U.b = as.matrix(vegdist(bioV[22:42,-1:-2], method = "bray"))
dbRDA.mat_Ub = capscale(bio_U.b ~ ., benthos_mt, comm = bioV[22:42,-1:-2])

# occurence
rich_e.jac = as.matrix(vegdist(rich_e, method = "jaccard"))
dbRDA.mat_eJ = capscale(rich_e.jac ~ ., benthos_mt, comm = rich_e)

vec_D = rep(0, 1000)
vec_U = rep(0, 1000)
vec_e = rep(0, 1000)

for (n in 1:1000) {
  # DOV
  df_D = sample(rich_D, replace = T) # bootstrap sampling
  sim_D = as.matrix(vegdist(df_D, method = "jaccard")) # similarity matrix
  mod_D = capscale(sim_D ~ ., benthos_mt, comm = df_D) # dbRDA
  vec_D[n] = RsquareAdj(mod_D)$r.squared # total constrained variation
  # UVC
  df_U = sample(rich_U, replace = T)
  sim_U = as.matrix(vegdist(df_U, method = "jaccard"))
  mod_U = capscale(sim_U ~ ., benthos_mt, comm = df_U)
  vec_U[n] = RsquareAdj(mod_U)$r.squared
  # eDNA
  df_e = sample(rich_e, replace = T)
  sim_e = as.matrix(vegdist(df_e, method = "jaccard"))
  mod_e = capscale(sim_e ~ ., benthos_mt, comm = df_e)
  vec_e[n] = RsquareAdj(mod_e)$r.squared
}
vec_D = sort(vec_D)
vec_U = sort(vec_U)
vec_e = sort(vec_e)

# The 25th and 975th of the bootstrapped values
vec_D[c(25, 975)]
vec_U[c(25, 975)]
vec_e[c(25, 975)]

tab = data.frame(RDA1_2 = vector(length = 7), 
                constrained = vector(length = 7),
                row.names = c("DOV_occur", "DOV_abun", "DOV_bio",
                              "UVC_occur", "UVC_abun", "UVC_bio", "eDNA"))
mod = list(dbRDA.mat_DJ, dbRDA.mat_DB, dbRDA.mat_Db,
           dbRDA.mat_UJ, dbRDA.mat_UB, dbRDA.mat_Ub,
           dbRDA.mat_eJ)
rda.plot = list()
col3 = c(rep("#F8766D", 3), rep("#619CFF", 3), "#00BA38")
shp3 = c(rep(c(16, 17, 15), 2), 16)
for (nu in 1:7) {
  smry = summary(mod[[nu]])
  tab$RDA1_2[nu] = round(sum(smry$cont$importance[2, 1:2])*100, 2) # constrained variations on RDA1&2
  tab$constrained[nu] = round(RsquareAdj(mod[[nu]])$r.squared, 2) # total constrianed variations
  df1  = data.frame(smry$sites[,1:2]) # site scores for RDA1 and RDA2
  df2  = data.frame(smry$biplot[,1:2])  # mapping benthic variables
  rda.plot[[nu]] = ggplot(df1, aes(x=CAP1, y=CAP2)) +
    ggtitle(paste(row.names(tab)[nu], "(", round(smry$cont$importance[3, 8]*100, 2), "% constrained)"))+
    geom_segment(data=df2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
                 color="grey50", arrow=arrow(length=unit(0.01,"npc"))) +
    geom_text(data=df2, aes(x=CAP1,y=CAP2,label=rownames(df2),
                            hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
              color="grey50", size=3) +
    geom_point(size = 3, col = col3[nu], shape = shp3[nu]) +
    geom_text(aes(label=rownames(df1),
                  hjust=0,vjust=1.5), colour = "black",size=3) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    xlab( paste("RDA1 (", round(smry$cont$importance[2, 1]*100, 2), "%)")) +
    # variations constrained in RDA1 
    ylab(paste("RDA2 (", round(smry$cont$importance[2, 2]*100, 2), "%)")) +
    # variations constrained in RDA2
    theme(plot.title = element_text(face = "bold", size = (15), hjust = 0.5),           axis.title = element_text(face = "bold", size = 14, colour = "black"),
          axis.text = element_text(colour = "black", face = "bold", size = 12), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "bottom", 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank())
}
# Figure 6. dbRDA on fish assemblages constrained by benthic composition.
ggarrange(rda.plot[[1]], rda.plot[[2]], rda.plot[[3]], rda.plot[[4]],
          rda.plot[[5]], rda.plot[[6]], rda.plot[[7]],
          nrow = 3, ncol = 3, legend = "bottom",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"))

# Extract genus from species name
spp = list(colnames(rich_D), colnames(rich_U), colnames(rich_e))
genus = vector(mode = "list", length = 3)
for (i in 1:3) {
  ge_sp = str_split(spp[[i]], pattern = " ")
  for (j in 1:length(ge_sp)) {
    genus[[i]][j] = ge_sp[[j]][1]
    genus[[i]] = row.names(table(genus[[i]]))}
  ge_sp = NULL
}

# Retreive taxonomic information from Fishbase
FB = rfishbase::load_taxa()
R3 = FB %>% filter(Genus %in% row.names(table(unlist(genus))))
RD = FB %>% filter(Genus %in% genus[[1]])
RU = FB %>% filter(Genus %in% genus[[2]])
Re = FB %>% filter(Genus %in% genus[[3]])

length(table(R3$Order)) # 36 orders
length(table(R3$Family)) # 106 families
length(table(R3$Genus)) # 241 genera
ncol(rich_3[-1:-2]) # 438 species

length(table(RD$Order)) # 10 orders
length(table(RD$Family)) # 18 families
length(table(RD$Genus)) # 46 genera
ncol(rich_D) # 83 species

length(table(RU$Order)) # 14 orders
length(table(RU$Family)) # 27 families
length(table(RU$Genus)) # 60 genera
ncol(rich_U) # 111 species

length(table(Re$Order)) # 36 orders
length(table(Re$Family)) # 106 families
length(table(Re$Genus)) # 233 genera
ncol(rich_e) # 383 species

setdiff(union(colnames(rich_D), colnames(rich_U)), colnames(rich_e))

bv = bc_mt[,-11]
avg = aggregate(x = bv/rowSums(bv)*100, by = list(NE$Area), FUN = function(x){round(mean(x),2)})
tb_avg = as.data.frame(t(avg[,-1]))
sdd = aggregate(x = bv/rowSums(bv)*100, by = list(NE$Area), FUN = function(x){round(sd(x),2)})
tb_sd = as.data.frame(t(sdd[,-1]))

colnames(tb_avg) = c("Northern coast", "Outlying islands")
colnames(tb_sd) = c("Northern coast", "Outlying islands")

tb = as.data.frame(matrix(nrow = 10, ncol = 4), row.names = colnames(bv))
colnames(tb) = c("Northern coast", "Outlying islands", "t-statistic", "p-value")

for (n in 1:10) {
  mod = t.test(bv[,n]~NE$Area)
  tb$`Northern coast`[n]= paste(tb_avg[n,1], tb_sd[n, 1], sep = " ± ")
  tb$`Outlying islands`[n]= paste(tb_avg[n,2], tb_sd[n, 2], sep = " ± ")
  tb$`t-statistic`[n] = round(mod$statistic, 2)
  tb$`p-value`[n] = round(mod$p.value, 3)
  }
tb

df_richness = as.data.frame(matrix(nrow = 21, ncol = 9),
                            row.names = NE$Code)
colnames(df_richness) = str_c(rep(c("DOV", "UVC", "eDNA"), each = 3),
                              rep(c("Species", "Genus", "Family"),3), sep = "_")
m3 = list(rich_D, rich_U, rich_e)
R3 = list(RD, RU, Re)
for (n in 1:3) {
  df = m3[[n]]
  df$Site = NE$Code
  df = melt(df, variable_name = "Species")
  df = left_join(df, R3[[n]][,c(2,3,5,6)]) # combine taxonomic information
  df = df[df$value == 1,] # remove species absent
  # species
  df_richness[,(1+3*(n-1))] = rowSums(m3[[n]])
  # genus
  df_richness[,(2+3*(n-1))] = df[,c(1,4)] %>% distinct() %>% group_by(Site) %>%
    summarise(Genus = n()) %>% pull(Genus)
  # Family
  df_richness[,(3+3*(n-1))] = df[,c(1,5)] %>% distinct() %>% group_by(Site) %>%
    summarise(Family = n()) %>% pull(Family)
}
# exporting Table
df_richness

df = data.frame(Species = as.numeric(unlist(df_richness[,c(1,4,7)])),
                Genus = as.numeric(unlist(df_richness[,c(2,5,8)])),
                Family = as.numeric(unlist(df_richness[,c(3,6,9)])),
                Method = as.factor(rep(c("DOV", "UVC", "eDNA"), each = 21)),
                Site = NE$Code)
df[,-5] %>% group_by(Method) %>% summarise_all(lst(mean, sd)) # mean & sd
# Species
friedman.test(y = df$Species, groups = df$Method, blocks = df$Site)
frdAllPairsConoverTest(y = df$Species, groups = df$Method, blocks = df$Site)
# Genus
friedman.test(y = df$Genus, groups = df$Method, blocks = df$Site)
frdAllPairsConoverTest(y = df$Genus, groups = df$Method, blocks = df$Site)
# Family
friedman.test(y = df$Family, groups = df$Method, blocks = df$Site)
frdAllPairsConoverTest(y = df$Family, groups = df$Method, blocks = df$Site)

rdf = as.data.frame(matrix(c(rowSums(rich_D), rowSums(rich_U), rowSums(rich_e)),
             nrow = 3, ncol = 21, byrow = T)) # species richness of all methods
row.names(rdf) = c("DOV", "UVC", "eDNA")

# plotting
col3 = c("#F8766D", "#619CFF", "#00BA38") # colors for DOV, UVC, and eDNA
gg = vector(mode = "list", length = 3)
vp = c(100, 100, 5)
for (n in 1:3) {
  df = data.frame(x = flow$Speed, y = as.numeric(rdf[n,]))
  mod = summary(lm(y ~ x, df))
  gg[[n]] = ggplot(df, aes(x, y)) +
    theme(plot.title = element_text(face = "bold", size = (15), hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_text(colour = "black", face = "bold", size = 12), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank())+
    ggtitle(row.names(rdf)[n])+
    geom_point(shape = 16, size = 3, col = col3[n]) +
    geom_smooth(formula = y~x, method = "lm")+
    geom_label(x = 0.125, y = vp[n], 
               label = paste("y=", round(mod$coefficients[2,1], 2),
                             "x+", round(mod$coefficients[1,1], 2), ", R2=",
                             round(mod$r.squared,2), ", p=", round(mod$coefficients[2,4], 2)))+
    scale_y_continuous(limits = c(0, 110))
}

ag = ggarrange(gg[[1]], gg[[2]], gg[[3]], nrow = 1, ncol = 3,
              labels = c("(a)", "(b)", "(c)"))
annotate_figure(ag, left = text_grob("Species richness", face = "bold", size = 18, rot = 90),
                bottom = text_grob("Flow speed (m/s)", face = "bold", size = 18))

abun = as.data.frame(cast(rbind(sp_D, sp_U), Method + Site ~ Species,
                          value='Number', fun.aggregate = sum))
gg = vector(mode = "list", length = 21)
for (i in 1:21) {
  R = data.frame(DOV = t(abun[i, -1:-2]),
                 UVC = t(abun[(i+21), -1:-2]))
  colnames(R) = c("DOV", "UVC")
  m = c(1, 10, 50, 100, 500, 1000, 2000)
  out = iNEXT(R, q = 0, datatype = "abundance", size = m)
  names(out$iNextEst) = colnames(R)
  gg[[i]] = ggiNEXT(out, type=1)+ggtitle(NE$Code[i])+
    scale_x_continuous(name = "Fish individuals") +
    scale_y_continuous(name = "Species richness") +
    scale_shape_manual(values= rep(20, 21))+
    theme(plot.title = element_text(face = "bold", size = (15), hjust = 0.5),
          axis.text = element_text(colour = "black", face = "bold", size = 12), 
          axis.title = element_blank(), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.position = "none", legend.key=element_blank())+
    labs(x = "Sampling units")}
ggr = ggarrange(gg[[1]], gg[[2]],  gg[[3]], gg[[4]], gg[[5]], gg[[6]], gg[[7]],
                gg[[8]], gg[[9]],  gg[[10]], gg[[11]], gg[[12]], gg[[13]], gg[[14]],
                gg[[15]], gg[[16]],  gg[[17]], gg[[18]], gg[[19]], gg[[20]], gg[[21]],
                nrow = 7, ncol = 3, common.legend = T)
annotate_figure(ggr, bottom = text_grob("Fish individuals", 
                                        color = "black", face = "bold", size = 16),
                left = text_grob("Species richness", 
                                 color = "black", face = "bold", size = 16, rot = 90))
