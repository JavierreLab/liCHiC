# Chromatin Assortativity (ChAs)

"This global measure identifies to what extent a property of a chromatin fragment is shared by fragments that interact preferentially with it" Pancaldi, V., Carrillo-de-Santa-Pau, E., Javierre, B.M. et al. Integrating epigenomic data and 3D genomic structure with a new measure of chromatin assortativity. Genome Biol 17, 152 (2016). [https://doi.org/10.1186/s13059-016-1003-3](https://doi.org/10.1186/s13059-016-1003-3)

With this kind of analysis we want to identify potential features (proteins, histone marks...) that could be somehow mediating the interactions.

## Summary of the workflow

1.  **Preparing Data**
2.  **Computing ChAs**: [chaser](https://bitbucket.org/eraineri/ChAseR/src)
3.  **Randomizations**: chaser


## 1. Preparing Data

The first step of the analysis is to prepare the data in the proper format to be analysed.

Basically, we need 2 files, one with the **interactions** to build our network, and another one with the feature to interrogate. For the interactions we use the ibed obtained from liCHi-C, and for the **feature** we use a bed file with the peak coordinates from ChIP-seq.

Due to the nature of our liCHi-C data we can explore 2 main subnetworks the Promoter-Promoter and the Promoter-OtherEnd (non-promoter regions). To deal with liCHi-C we developed an in-house R package [HiCaptuRe](https://github.com/LaureTomas/HiCaptuRe).

## 2. Computing ChAs

```R

## we define 2 variables, for the interactions file and for the feature file(s)
ibed <- "path/to/interactions.ibed"
features <- c("path/to/H3K4me3.bed","path/to/H3K4me1.bed")
names(features) <- c("H3K4me3","H3K4me1")

# we load the ibed file using the load_interactions function from HiCapture
interactions <- HiCaptuRe::load_interactions(ibed)

# we retrieve the coordinates of both ends for the 3 networks: complete, and P-P and P-OE subnetworks
complete_network <- as.data.frame(interactions)[,c(1,2,3,10,11,12)]
subnetwork_PP <- as.data.frame(interactions[interactions$int == "P_P"])[,c(1,2,3,10,11,12)]
subnetwork_POE <- as.data.frame(interactions[interactions$int == "P_OE"])[,c(1,2,3,10,11,12)]

# then using chaser we make the network objects for each of them
net <- chaser::make_chromnet(complete_network)
netPP <- chaser::make_chromnet(subnetwork_PP)
netPOE <- chaser::make_chromnet(subnetwork_POE)
  
# we load the feature(s) to the networks
for (number in 1:length(features)) 
{
  net <- chaser::load_features(net, features[number], type="bed3",
                               missingv=0,featname = names(features)[number])
  netPP <- chaser::load_features(netPP, features[number], type="bed3",
                                 missingv=0,featname = names(features)[number])
  netPOE <- chaser::load_features(netPOE, features[number], type="bed3",
                                 missingv=0,featname = names(features)[number])
}

# and finally we compute the ChAs value for all the features in each network at once
net_chas <- chaser::chas(net)
net_chasPP <- chaser::chas(netPP)
net_chasPOE <- chaser::chas(netPOE)

# Additionally we can compute the abundancy for all features in each network, to use for visualization for example
abund <- colMeans(chaser::export(net = net))
abundPP <- colMeans(chaser::export(net = netPP))
abundPOE <- colMeans(chaser::export(net = netPOE))
```

## 3. Randomization

This step is esential to validate if our ChAs value is significant or not. Inside the Chaser R package there are 3 types of randomizations: vainilla, preserve.nodes and dist.match. In our case we use the dist.match to maintain the distance of the interactions, since it's a specific characteristic of our networks.

```R
# we use the randomize function with at least 1000 randomizations
rnets <- chaser::randomize(net, nrandom=1000, dist.match=TRUE)
rnetsPP <- chaser::randomize(netPP, nrandom=1000, dist.match=TRUE)
rnetsPOE <- chaser::randomize(netPOE, nrandom=1000, dist.match=TRUE)
```

Finally once we have out ChAs values and the associate randomizations values we can test if our value is significant.
