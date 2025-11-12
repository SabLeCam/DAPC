# DAPC

support de cours sur la DAPC, basé sur un cours d'Amélia Viricel

Documentation source:

Jombart, T., Devillard, S. & Balloux, F. Discriminant analysis of principal components: a new method for the analysis of genetically structured populations. BMC Genet 11, 94 (2010). [https://doi.org/10.1186/1471-2156-11-94](https://link.springer.com/article/10.1186/1471-2156-11-94#article-info)

Cullingham, C., Peery, R. M., & Miller, J. M. (2023). A roadmap to robust discriminant analysis of principal components. Molecular ecology resources, 23(3), 519-522. [https://doi.org/10.1111/1755-0998.13724](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13724?casa_token=_nnhWtUcXrsAAAAA%3AnsTgaLQIcI1tRiKRD1pkxSDwQxP9cDvH5E9Yof78xrfLWWQoGm3NLYauqs2dSOxI08vTfUCfMAGUuZOS)

## Qu'est ce que la DAPC (*Discriminant Analysis of Principal Components*)?

Une méthode multivariée conçue pour identifier et décrire des groupes d'individus génétiquement apparentés. Elle combine les principes de l'ACP (Analyse en Composante Principale) et l'Analyse Discriminante.
L'ACP vise à synthétiser la variabilité globale entre les individus, incluant à la fois la divergence entre les groupes (variabilité génétique structurée) et la variation intragroupe (variabilité génétique aléatoire). Pour évaluer les relations entre différents groupes, une méthode adéquate doit se concentrer sur la variabilité intergroupe, en négligeant la variation intragroupe. C'est précisément le principe de l'analyse discriminante (AD). Elle définit un modèle dans lequel la variation génétique est partitionnée en une composante intergroupe et une composante intragroupe, et produit des variables synthétiques qui maximisent la première tout en minimisant la seconde. En d'autres termes, l'analyse discriminante (AD) tente de résumer la différenciation génétique entre les groupes, tout en négligeant la variation intragroupe.

Figure issue de Jombart et al. 2010 <img width="445" height="861" alt="image" src="https://github.com/user-attachments/assets/28a7767d-8816-42c9-9b48-d6faf6dc171e" />

## Comment utiliser la DAPC?

Cette méthode d'analyse a l'avantage d'être beaucoup moins demandante en ressources computationnelles que les méthodes de clustering bayésiennes, comme STRUCTURE. Elle est donc communément utilisée pour réaliser une première analyses de la distribution de la variabilité génétique dans le jeu de données. Dans ce cas , elle peut être utilisée __*sans à priori*__. c'est à dire qu'on ne fait pas d'hypothèse quant au nombre de cluster dans notre jeu de données mais on se base sur des paramètres statistiques (BIC,Bayesian Information Criterion issus de comparaisons de modèles) pour estimer le nombre le plus probable de groupes génétique.
Ainsi dans une première étape, on estime le nombre de clusters génétique puis dans un deuxième temps, on estime la probablilité d'assignation des individus à ces différents clusters.

Elle peut aussi être utilisée avec __*à priori*__. Dans ce cas, des informations sont connues quand au nombre de populations probables dans le jeu de données, et cela permet de valider ou non cette hypothèse. C'est aussi utile dans le cas où la migration entre les clusters est forte. Appliquer une nombre de cluster de  K (nb pop) - 1 permet une meilleur estimation de l'assignation des individus aux différents clusters.

## PREMIER EXEMPLE: *Delphinus delphis*
https://whalesense.org/wp-content/uploads/2020/07/Dd-breach-1-1024x638.jpg<img width="512" height="319" alt="image" src="https://github.com/user-attachments/assets/a192d180-e445-414a-bdf3-2e243319a648" />

Jeu de données d' Amélia Viricel: 
83 individus échantillonnés en Altantique Nord Est (Golfe de Gascogne) ont été séquencés avec les méthode du sdRADseq (NotI)
Après filtration, 10, 930 SNP ont été conservés.

Nous allons réaliser une DAPC sur ce jeu de données.

Nous allons utiliser le package *adegenet*
```R
#install.packages("adegene")  # si necessaire
library(adegenet)
```
Importer les données sous forme Genepop
```R
C <- read.genepop("RAD_dauphins_AV.gen", quiet = TRUE)
C 				## to make sure the data look OK (nb of individuals and loci)
```
Etape 1: estimer le nombre de cluster qui décrit le mieux la variablité génétique présente dans les données (ici, le nb de PC testés est 40)
```R
grp <- find.clusters(C, n.pca=40)
```
Cette fonction va créer une figure représentant les valeurs du BIC. De façon itérative, on décide combien de Clusters (K) on retient (la plus faible valeur de BIC)

```R
grp$size  		## will return how many individuals are assigned to each cluster
#write.table(head(grp$grp,374), file = "clusters.csv", sep = ",", col.names = NA)  ## to export membership to each cluster as csv file
```
Maintenant on réalise la DAPC avec le nb de clusters choisis
Il faut d'abord choisir le nb de PC retenue. Il faut en retenir suffisamment pour 80% de la variance soit expliquée. Ici 50
Ensuite on choisi le nb de fonctions disciminantes. Lorsque que le nb de clusters est petit, on peut toutes les garder.
```R
dapc1 <- dapc(C, grp$grp)   
dapc1
#write.table(dapc1$posterior, file="dapc1_MembershipProb", sep=",", row.names=T) ## pour faire une table des membership Prob
```
On peut alors représenter les probablités d'assignation des individus aux différents clusters.
```R
compoplot(dapc1, posi="bottomright", txt.leg=paste("Cluster", 1:2), lab="", ncol=2, xlab="individuals") ## pour plotter les membership proba
```
<img width="685" height="467" alt="image" src="https://github.com/user-attachments/assets/9efbed57-d851-4e62-86aa-ec12cc3f454d" />

- Que pouvez-vous déduire de ces analyses?
- Quel a été l'intérêt de la DAPC à ce stade de l'analyse?

## DEUXIEME EXEMPLE: *Salmo salar*
Ce jeu de données est issu de l'article suivant:

Moore, J. S., Bourret, V., Dionne, M., Bradbury, I., O'Reilly, P., Kent, M., ... & Bernatchez, L. (2014). Conservation genomics of anadromous Atlantic salmon across its North American range: outlier loci identify the same patterns of population structure as neutral loci. Molecular Ecology, 23(23), 5680-5697. *Disponible sur le Moodle du cours.*

Carte d'échantillonnage:
<p float="left">
<img width="284" height="398" alt="image" src="https://github.com/user-attachments/assets/6b4fcc84-6e97-4c01-9c45-6ac250f66081" />
 <img width="400" height="160" alt="image" src="https://github.com/user-attachments/assets/e262476e-2034-4a4c-8cfd-20d47fe86313" />
</p>
Sources: Moore et al. 2010, © USFWS.

Le jeu de données est constitué de 1,080 échantillons génotypes à 3092 SNP à l'aide d'une puce à SNP. Les échantillons ont été prélevés dans 49 sites. 

Importer les données

```R
C <- read.genepop(ncode = 3L, "Ssalar_SNPall_genepop_28Apr2014.gen", quiet = TRUE) ## to import your data from a Genepop file
C 				## to make sure the data look OK (nb of individuals and loci)
summary(C)      ## to look at some info about your data (nb of alleles, missing data, etc)
```
Estimer le nb de clusters
```R
grp <- find.clusters(C, max.n.clust=40)   ## function to find the number of clusters that best describe the variability in the data (here with a max nb of clusters tested of 40)
## function will plot percent variance explained vs nb of PCs and will ask how many PCs should be retained: retain all PCs at this step, then looking at the plot of BIC values, decide how many clusters are best (lowest BIC) and enter K chosen (K=9)
grp$size  		## will return how many individuals are assigned to each cluster :  [1]  99 144  68  39 136 231 170  75  47  71
write.table(head(grp$grp,374), file = "clusters.csv", sep = ",", col.names = NA)  ## to export membership to each cluster as csv file; 1080 is the number of inviduals
```
On a choisi K=9. On reéalise donc l'analyse avec K = 9
```R
dapc1 <- dapc(C, grp$grp)    ## to perform DAPC, function will ask how many PCs to retain : chose nb of PCs to retain about 80% for variance, then function asks how many discriminant functions should be examined, can keep all for small number of clusters (here 5)
dapc1   ### will return info about DAPC analysis
```R
scatter(dapc1, posi.da="topright")			## will do scatterplot of all clusters
```
<img width="685" height="467" alt="image" src="https://github.com/user-attachments/assets/ce2b7e3b-c9bc-4600-be0b-fc3f74e9acc3" />

```R
#write.table(dapc1$posterior, file="dapc1_MembershipProb", sep=",", row.names=T) ## pour faire une table des membership Prob
compoplot(dapc1, posi="bottomright", txt.leg=paste("Cluster", 1:9), lab="", ncol=2, xlab="individuals") ## pour plotter les membership proba
```
<img width="685" height="467" alt="image" src="https://github.com/user-attachments/assets/3ca1813e-86da-46de-b527-5ab8c74b6cbc" />


