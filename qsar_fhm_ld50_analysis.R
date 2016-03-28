#qsar_fhm_ld50_analysis.R
#GHS Chemicals: http://www.unece.org/fileadmin/DAM/trans/danger/publi/ghs/ghs_rev05/English/04e_part4.pdf
#See table 4.1.1
#Had to convert smiles to SDF on command line using: obabel fhm_chems2.smi -O fhm_chems2.sdf --gen2d

library(ChemmineR)
library(ChemmineOB)
#library(ggfortify)
#library(cluster)
#library(caret)
#library(randomForest)
library(arules)

distance_function <- function(x){
  fpset <- as(x, "FPset")
  simMAap <- sapply(cid(fpset), function(x) fpSim(x=fpset[x], fpset, sorted=FALSE))
  return(as.dist(1-simMAap))
}

decipher_ap_code <- function(ap_set, ap_code, compound){
  chem_row <- which(cid(apset) == compound)
  return(db.explain(apset[chem_row])[which(ap(apset[[chem_row]]) == ap_code)])
}

fhm_data <- read.delim(file="fhm_chems.sml", sep="\t", header=TRUE)
#fhm_category2 <- rep("Toxic", length(fhm_data$FHM_LD50))
#fhm_category2[which(fhm_data$FHM_LD50 > 50)] <- "NonToxic"

#QSAR Atom-Pair Fingerprinting
sdfset <- read.SDFset("fhm_chems2.sdf")
apset <- sdf2ap(sdfset)
fpset <- desc2fp(apset, descnames=4096, type="FPset")
simMAap <- sapply(cid(fpset), function(x) fpSim(x=fpset[x], fpset, sorted=FALSE))
rownames(simMAap) <- fhm_data$Name
colnames(simMAap) <- fhm_data$Name
hc <- hclust(as.dist(1-simMAap), method="single")
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=FALSE)
plot(hc)

#PCA
pc.cr <- prcomp(t(as.matrix(fpset)))
plot(pc.cr$rotation[,1], pc.cr$rotation[,2])
text(pc.cr$rotation[,1], pc.cr$rotation[,2], labels = fhm_data$Name)

#Fanny clustering
autoplot(fanny(t(as.matrix(fpset)), 4, memb.exp = 1.2), frame = TRUE)

#QSAR FHM LD50
hazard_cat <- rep("LessHazardous", nrow(fhm_data))
quantile(fhm_data$FHM_LD50, probs=c(0.4, 0.5, 0.8))
hazard_cat[which(fhm_data$FHM_LD50 <= 50)] <- "Hazardous"
#filter out columns that have no entries
col_sums <- colSums(fpset@fpma)
fpma_filtered <- fpset@fpma[, which(col_sums > 0)]


obj <- tune(svm, train.x=fpma_filtered, train.y=as.factor(hazard_cat), 
            ranges = list(gamma = 2^(-1:1), cost = 2^(2:4)),
            tunecontrol = tune.control(sampling = "fix")
)

svm_model <- svm(x=fpma_filtered, y=as.factor(hazard_cat), kernel="radial", type="C-classification", cost=100)
pred <- predict(svm_model, fpma_filtered)
table(pred, as.factor(hazard_cat))

#RandomForest
fhm.rf <- randomForest(x=fpma_filtered, y=as.factor(hazard_cat), importance=TRUE, proximity = TRUE, ntree = 2000, mtry = 9)
print(fhm.rf)
fhm.rf.regression <- randomForest(x=fpma_filtered, y=fhm_data$FHM_LD50, importance=TRUE, proximity = TRUE, mtry=10, ntree=2000)
print(fhm.rf.regression)

#Apriori algorithm
#As expected, this gives some better results
#And it allows for a Bayesian approach on the back-end
fpma_tox_levels <- cbind(fpma_filtered, CAT1=fhm_data$CAT1, CAT2=fhm_data$CAT2, CAT3=fhm_data$CAT3, CAT4=fhm_data$CAT4)
rules <- apriori(fpma_tox_levels, parameter = list(conf=0.50, target = "rules"))

rules.sub1 <- subset(rules, subset = rhs %pin% "CAT1")
inspect(sort(rules.sub1))
rules.sub2 <- subset(rules, subset = rhs %pin% "CAT2")
inspect(sort(rules.sub2, by="confidence"))
rules.sub3 <- subset(rules, subset = rhs %pin% "CAT3")
inspect(sort(rules.sub3, by="confidence"))
rules.sub4 <- subset(rules, subset = rhs %pin% "CAT4")
inspect(sort(rules.sub4, by="confidence"))

lhs_rule2.sub <- as(subset(rules, subset = lhs %pin% "54898248832" & lhs %pin% "54897200256" & lhs %pin% "54896151680" & rhs %in% "CAT2")@lhs, "matrix")

itemFrequencyPlot(as(fpma_tox_levels, "transactions"))


decipher_ap_code(ap_set=apset, ap_code="69933794434", compound = "CMP2")

#Cross-validation
#rules <- apriori(fpma_tox_levels, parameter = list(conf=.75, target = "rules"))

chem_fingerprints_matrix <- fpma_filtered[, which(colSums(fpma_filtered) > 0)]

fpma_filtered <- fpset@fpma[, which(col_sums > 0)]


#x is fpma_tox_levels
#ground_truth is a vector of the ground truth labels in same order as chems in x or fhm_data$GHS.Category in this case
apriori_cross_validation <- function(x, ground_truth, chem_fingerprints_matrix, folds = 10, conf=0.75, target="rules"){
  category_levels <- levels(ground_truth)
  chem_final_category <- NULL
  cross_validation_folds <- createFolds(fhm_data$GHS.Category, folds)
  for(i in 1:length(cross_validation_folds)){
    print(paste("i = ", i))
    rules <- apriori(x[-cross_validation_folds[[i]],], parameter = list(conf=.75, target = "rules"))
    sub.rules <- subset(rules, subset = rhs %pin% "CAT")
    rules.df <- as(sub.rules, "data.frame")
    chem_category_predicted <- NULL
    sub.rules <- NULL
    priors <- as.data.frame(table(ground_truth))
    prior_prob <- priors$Freq / length(ground_truth)
    priors <- cbind(priors, prior_prob = prior_prob)
    #for the kth chemical...
    for(k in 1:nrow(x)){
      chem_fingerprint <- names(chem_fingerprints_matrix[k,])
      chem_rules_stats <- NULL
      for(j in 1:length(category_levels)){
        print(paste("j = ", j))
        #Do this for each CAT, so CAT1, CAT2, CAT3, CAT4, CAT5
        sub.rules <- subset(rules, subset = rhs %pin% category_levels[j])
        rule_fits_compound <- apply(as(sub.rules@lhs, "matrix"), 1, FUN=rule_filter_in_compound, chem_fingerprint)
        whole_rules_fits_compound <- as(sub.rules, "data.frame")[which(rule_fits_compound == TRUE),]
        best_predictions <- whole_rules_fits_compound[which(whole_rules_fits_compound$confidence == max(whole_rules_fits_compound$confidence)),]
        best_predictions_empty <- nrow(best_predictions)
        if(best_predictions_empty > 0){
          posterior_prob <- best_predictions$confidence * priors$prior_prob[which(priors$ground_truth == category_levels[j])]
          chem_rules_stats_temp <- cbind(run = i, chemical_name = rownames(x)[k], category = category_levels[j], best_predictions, posterior_probability = posterior_prob)
          chem_rules_stats <- rbind(chem_rules_stats, chem_rules_stats_temp)
        }
      }
      chem_final_category_temp <- chem_rules_stats[which(chem_rules_stats$posterior_probability == max(chem_rules_stats$posterior_probability)),]
      chem_final_category <- rbind(chem_final_category, chem_final_category_temp)
    }
  }
  return(chem_final_category)
}
start_time <- proc.time()
x <- apriori_cross_validation(fpma_tox_levels, fhm_data$GHS.Category, chem_fingerprints_matrix, conf=0.75, target="rules")
proc.time() - start_time

find_max <- function(x){
  if(which(x$category == "CAT2") > 0){
    print("Yay")
  }
  x_max_post_prob <- first(x[which(x$posterior_probability == max(x$posterior_probability)),])
  return(x_max_post_prob)
}

xy <- by(x, x[,c("run", "chemical_name", "category")], find_max)
xy[["1", "CMP2", "CAT3"]]


