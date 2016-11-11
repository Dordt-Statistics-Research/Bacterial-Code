root_dir = "/home/jason/Documents/School/College/Fall Junior/STAT RESEARCH/MyWD/Bacterial-Code/"
setwd(root_dir)

load("bayesian3_exploratory_inputs/eijs.Rdata")
load("bayesian3_exploratory_inputs/operons.Rdata")
load("bayesian3_exploratory_inputs/uniaijs.Rdata")  # aijs from the UniMM-MI method
load("bayesian3_exploratory_inputs/aijResults.Rdata")  # Includes both aijs and fits from the MultiMM-MI method
multiaijs <- aijResults$aijs
source("expnameMaps.R")
source("Aij_utility_funcs.R")

cel_chip_map = get.expname.map("inputs/Ecoli_Cel_to_chip_fromClaire_edited.tab")

exp_data_for <- function(gene) as.matrix(eijs)[gene,]

hist_gene <- function(gene, ...) hist(exp_data_for(gene), breaks=25, ...)

get_operon <- function(gene) {
  which(sapply(operons, function(op) gene %in% op))
}

gene <- function(pegnum) paste0("fig|83333.1.peg.",pegnum)

dnormOff <- function(gene, x) (1-aijResults$pi[gene])*dnorm(x, mean=aijResults$mu0[gene], sd=sqrt(aijResults$var[[get_operon(gene)]][gene,gene]))
dnormOn <- function(gene, x) aijResults$pi[gene]*dnorm(x, mean=aijResults$mu1[gene], sd=sqrt(aijResults$var[[get_operon(gene)]][gene,gene]))

hist_with_overlay <- function(pegnum,figure) {
  hist_gene(gene(pegnum), main=paste(figure,"Expression data and predicted activity"), xlab="Log expression value", ylab="Density", prob=T)
  curve(dnormOff(gene(pegnum),x), add=T, col="darkblue", lwd=3)
  curve(dnormOn(gene(pegnum),x), add=T, col="darkred", lwd=3)
  legend("topright", c("Inactive prediction", "Active prediction"), fill=c("darkblue", "darkred"))
}

hist_with_overlay_gene <- function(gene,figure,range,...) {
  fig_id <- com_to_fig(gene)
  off <- list(seq(0,16,0.1),(dnormOff(fig_id,seq(0,16,0.1))))
  on <- list(seq(0,16,0.1),(dnormOn(fig_id,seq(0,16,0.1))))
  names(off) <- c("x","y")
  names(on) <- c("x","y")
  plot(off$x,off$y, col="darkblue", lwd=5,type="l",main="", xlab="", ylab=paste(gene,"Density"), xlim=range,bty="n",...)
  lines(on$x,on$y, col="darkred", lwd=5)
}

hist_with_overlay_gene_vert <- function(gene,figure,range,...) {
  fig_id <- com_to_fig(gene)
  off <- list(seq(0,16,0.1),(dnormOff(fig_id,seq(0,16,0.1))))
  on <- list(seq(0,16,0.1),(dnormOn(fig_id,seq(0,16,0.1))))
  names(off) <- c("x","y")
  names(on) <- c("x","y")
  plot(off$y,off$x, col="darkblue", lwd=5,type="l",main="", xlab=paste(gene,"Density"), ylab="", ylim=range,bty="n",...)
  lines(on$y,on$x, col="darkred", lwd=5)
}

#Save histograms overlays for all transcription factors in regprecise
save_hists <- function() {
  trials = c(2761, 4001, 310, 1311, 3292, 3357, 673, 4298, 826, 3755, 3414, 3703, 3788, 1173, 483, 1644, 2669, 2501, 81, 2664, 163, 4088, 113, 3956, 3883, 410, 3927, 342, 3170, 747, 2643, 3691, 4151, 1606, 3370, 1001, 3540, 4301, 460, 502, 1517, 803, 3973, 3953, 3181, 3889, 1811, 666, 3022, 20, 3226, 3508, 3076, 1837, 1328, 3041, 4234, 2167, 4014, 1323, 2764, 1389, 2792, 2372, 3826, 3827, 77, 2080, 4101, 2662, 3857, 2395, 1499, 65, 1308, 3792, 1528, 4251, 1425)
  #trials = c(2761, 4001, 310, 1311, 3292, 3357, 673, 4298, 826, 3755, 3414, 3703, 3788, 1173, 483, 1644, 2669, 2501, 81, 2664, 163, 4088, 113, 3956, 3883, 410, 3927, 342, 3170, 747, 2643, 3691, 4151, 1606, 3370, 1001, 3540, 4301, 460, 502, 1517)
  
  for (trial in trials) {
    png(filename = paste0("/media/jason/PNY_BLUE/STATS/Meeting 10-17-16/TFs/", as.character(trial)))
    hist_with_overlay(trial, trial)
    dev.off()
  }
}

basicplot <- function(pegnum1, pegnum2,aijs,method,gene1,gene2,xlimit,ylimit) {
  plot(exp_data_for(gene(pegnum1)), exp_data_for(gene(pegnum2)),main=method,col="black",xlim=xlimit,ylim=ylimit,pch=16,xlab=gene1,ylab=gene2);
}

eij_vs_eij <- function(pegnum1, pegnum2,aijs,method,gene1,gene2,xlimit,ylimit,tf_and_target=FALSE,pegnum3=-1,gene3="") {
  eij_vs_eij_helper(pegnum1, pegnum2,aijs[gene(pegnum1),],aijs[gene(pegnum2),],method,gene1,gene2,xlimit,ylimit,tf_and_target,pegnum3,aijs[gene(pegnum3),],gene3)
}

eij_vs_eij_helper <- function(pegnum1,pegnum2,aijs1,aijs2,method,gene1,gene2,xlimit,ylimit,tf_and_target=FALSE,pegnum3=-1,aijs3=c(),gene3="") {
  
  xlimit = c(min(exp_data_for(gene(pegnum1)))-1, max(exp_data_for(gene(pegnum1)))+1)
  ylimit = c(min(exp_data_for(gene(pegnum2)))-1, max(exp_data_for(gene(pegnum2)))+1)
  
  inact_hi1=which(aijs1<0.2);
  act_hi1=which(aijs1>=0.8);
  other1=intersect(which(aijs1>=0.2),which(aijs1<0.8));
  
  inact_hi2=which(aijs2<0.2);
  act_hi2=which(aijs2>=0.8);
  other2=intersect(which(aijs2>=0.2),which(aijs2<0.8));
  
  inact_hi3=which(aijs3<0.2);
  act_hi3=which(aijs3>=0.8);
  other3=intersect(which(aijs3>=0.2),which(aijs3<0.8));
  
  conflicts=NULL;
  other=NULL;

  #if pegnum1/gene1 is a transcription factor and pegnum2/gene2 is a target of that transcription factor
  if (tf_and_target) {
    conflicts=intersect(act_hi1, inact_hi2)
    if (length(conflicts) > 0) {
      writeLines(paste0("Conflict experiments for ", gene1, " and ", gene2, " in method ",  method, ": "))
      cel_names = (names(eijs)[conflicts])
      for (i in c(1:length(cel_names))) {
        cel_names[i] = substr(cel_names[i], 1,nchar(cel_names[i])-7)
      }
      lapply(get_chipname(cel_chip_map,cel_names), FUN=function(x)(writeLines(paste0("\t",x))))
    }
    other=setdiff(which(is.finite(aijs1)), conflicts)
    #par(mfrow=c(2,3))
  }else{
    #par(mfrow=c(2,2))
  }
  
  sym_cex = 0.75
  
  plot(c(),c(),main=method,xlim=xlimit,ylim=ylimit,xlab=gene1,ylab=gene2);
  points(exp_data_for(gene(pegnum1))[act_hi1], exp_data_for(gene(pegnum2))[act_hi1],col="darkred",pch=17,cex=sym_cex)
  points(exp_data_for(gene(pegnum1))[inact_hi1], exp_data_for(gene(pegnum2))[inact_hi1],col="darkblue",pch=16,cex=sym_cex)
  points(exp_data_for(gene(pegnum1))[other1], exp_data_for(gene(pegnum2))[other1],col="black",pch=0,cex=sym_cex)
  legend("topleft", c(paste(gene1, "Active\t\t", length(act_hi1)), paste(gene1, "Inactive\t", length(inact_hi1)), paste("Unsure\t\t\t", length(other1))),pch=c(17,16,0),col=c("darkred","darkblue","black"),pt.cex=1,cex=0.8)
  
  plot(c(),c(),main=method,xlim=xlimit,ylim=ylimit,xlab=gene1,ylab=gene2);
  points(exp_data_for(gene(pegnum1))[act_hi2], exp_data_for(gene(pegnum2))[act_hi2],col="darkorange",pch=17,cex=sym_cex)
  points(exp_data_for(gene(pegnum1))[inact_hi2], exp_data_for(gene(pegnum2))[inact_hi2],col="darkgreen",pch=16,cex=sym_cex)
  points(exp_data_for(gene(pegnum1))[other2], exp_data_for(gene(pegnum2))[other2],col="black",pch=0,cex=sym_cex)
  legend("topleft", c(paste(gene2, "Active\t\t", length(act_hi2)), paste(gene2, "Inactive\t", length(inact_hi2)), paste("Unsure\t\t\t", length(other2))),pch=c(17,16,0),col=c("darkorange","darkgreen","black"),pt.cex=1,cex=0.8)
  
  if (!(pegnum3==-1)) {
    plot(c(),c(),main=method,xlim=xlimit,ylim=ylimit,xlab=gene1,ylab=gene2);
    points(exp_data_for(gene(pegnum1))[act_hi3], exp_data_for(gene(pegnum2))[act_hi3],col="darkslateblue",pch=17,cex=sym_cex)
    points(exp_data_for(gene(pegnum1))[inact_hi3], exp_data_for(gene(pegnum2))[inact_hi3],col="darkslategray2",pch=16,cex=sym_cex)
    points(exp_data_for(gene(pegnum1))[other3], exp_data_for(gene(pegnum2))[other3],col="black",pch=0,cex=sym_cex)
    legend("topleft", c(paste(gene3, "Active\t\t", length(act_hi3)), paste(gene3, "Inactive\t", length(inact_hi3)), paste("Unsure\t\t\t", length(other3))),pch=c(17,16,0),col=c("darkslateblue","darkslategray2","black"),pt.cex=1,cex=0.8)
    
  }
  
  if (tf_and_target) {
    plot(c(),c(),main=method,xlim=xlimit,ylim=ylimit,xlab=gene1,ylab=gene2);
    points(exp_data_for(gene(pegnum1))[conflicts], exp_data_for(gene(pegnum2))[conflicts],col="black",pch=13,cex=sym_cex)
    points(exp_data_for(gene(pegnum1))[other], exp_data_for(gene(pegnum2))[other],col="darkgrey",pch=20,cex=sym_cex)
    legend("topleft", c(paste("TRN Conflicts\t", length(conflicts)), paste("Other\t\t\t", length(other))),pch=c(13,20),col=c("black","darkgrey"),pt.cex=1,cex=0.8)
  }
}

trn <- function(tf,target,aijs,method,gene1,gene2,xlimit,ylimit,third=-1,gene3="") {
  
  on_threshold = 0.8 #genes with aijs greater than or equal to this will be considered to be on
  off_threshold = 0.2 #genes with aijs less than this will be considered to be off
  
  tf_aijs = aijs[gene(tf),]
  target_aijs = aijs[gene(target),]
  if (! third==-1) {
    third_aijs = aijs[gene(third),]
  }
  
  new_tf_aijs = tf_aijs
  new_target_aijs = target_aijs
  
  tf_active = which(tf_aijs>=on_threshold)
  target_inactive = which(target_aijs<off_threshold)
  
  new_target_aijs[tf_active] = 0.5*tf_aijs[tf_active] + 0.5*target_aijs[tf_active]
  new_tf_aijs[target_inactive] = 0.5*tf_aijs[target_inactive] + 0.5*target_aijs[target_inactive]
  
  eij_vs_eij_helper(tf,target,new_tf_aijs,new_target_aijs,method,gene1,gene2,xlimit,ylimit,TRUE,third,third_aijs,gene3)
}

compare_genes_helper <- function(tf, target, tf_name, target_name, third=-1, third_name="") {
  if (third==-1) {
    par(mfrow=c(2,3))
  } else {
    par(mfrow=c(2,4))
  }
  
  eij_vs_eij(tf,target,multiaijs,"MultiMM",tf_name,target_name,xlimit=c(2,16),ylimit=c(2,16),TRUE,third,third_name)
  trn(tf,target,multiaijs,"TRN_Multi",tf_name,target_name,xlimit=c(2,16),ylimit=c(2,16),third,third_name)
  
  #eij_vs_eij(tf,target,uniaijs,"UniMM",tf_name,target_name,xlimit=c(2,16),ylimit=c(2,16),TRUE,third,third_name)
  #trn(tf,target,uniaijs,"TRN_Uni",tf_name,target_name,xlimit=c(2,16),ylimit=c(2,16),third,third_name)
}

compare_genes <- function(map, tf_name, target_name, third_name) {
  a = get_pegnum(tf_name, map)
  b = get_pegnum(target_name, map)
  c = tf_name
  d = target_name
  e = get_pegnum(third_name, map)
  f = third_name
  compare_genes_helper(a, b, c, d, e, f)
}

get_pegnum <- function(gene_name, map) {
  #get the list of names from the gene map
  #get the peg id
  #split the ped id by the periods
  #flatten the list
  #return the fourth element, the actual numeric part
  s = unlist(strsplit(map[[gene_name]][["peg"]], split="[.]"))
  return(s[4])
}



plot.elaborate <- function(gene1,gene2,low,high) {
  fig1 <- com_to_fig(gene1)
  fig2 <- com_to_fig(gene2)
  
  data1 <- exp_data_for(fig1)
  data2 <- exp_data_for(fig2)
  
  xlimit = c(min(data1)-1, max(data1)+1)
  ylimit = c(min(data2)-1, max(data2)+1)
  
  aijs1 <- multiaijs[fig1,]
  aijs2 <- multiaijs[fig2,]
  
  inact1=which(aijs1<low);
  act1=which(aijs1>=high);
  other1=intersect(which(aijs1>=low),which(aijs1<high));
  
  inact2=which(aijs2<low);
  act2=which(aijs2>=high);
  other2=intersect(which(aijs2>=low),which(aijs2<high));
  
  regions <- list(intersect(inact1,inact2),
                  intersect(other1,inact2),
                  intersect(act1,inact2),
                  intersect(inact1,other2),
                  intersect(other1,other2),
                  intersect(act1,other2),
                  intersect(inact1,act2),
                  intersect(other1,act2),
                  intersect(act1,act2))
  
  colors <- c("darkblue","darkgrey","black","darkgrey","darkgrey","darkgrey","forestgreen","darkgrey","darkred")
  symbols <- c(16,18,4,18,18,18,15,18,17)
  
  layout(matrix(c(2,4,1,3),ncol=2,byrow=TRUE),heights=c(0.3,0.7), widths=c(0.7,0.3))
  
  sym_cex = 1.25
  font_cex = 1.75
  
  par(oma=c(0,0,2,0), mar=c(5,5,0,0))
  plot(c(),c(),main="",xlim=xlimit,ylim=ylimit,xlab=paste(gene1,"Expression Levels"),ylab=paste(gene2,"Expression Levels"),cex.lab=font_cex,cex.axis=font_cex)
  for (i in 1:length(regions)) {
    x_cex = sym_cex
    if (symbols[i]==4) {
      x_cex = sym_cex+1.25
    }
    points(data1[regions[[i]]], data2[regions[[i]]],col=colors[i],pch=symbols[i],cex=x_cex)
  }
  #box("plot", col="red")
  #box("figure", col="blue")
  
  
  par(mar=c(0,5,2,0)) 
  hist_with_overlay_gene(gene1,"",xlimit,cex.lab=font_cex,xaxt="n",yaxt="n",cex.axis=font_cex)
  #box("plot", col="red")
  #box("figure", col="blue")
  
  
  par(mar=c(5,0,0,2))
  hist_with_overlay_gene_vert(gene2,"",ylimit,cex.lab=font_cex,xaxt="n",yaxt="n",cex.axis=font_cex)
  #box("plot", col="red")
  #box("figure", col="blue")
  
  par(mar=c(0.1,0.1,2,0.1)) 
  plot.new()
  legend("center",c(paste(gene1,"inact.,",gene2,"inact."),
                  paste(gene1,"inact.,",gene2,"act."),
                  paste(gene1,"act.,",gene2,"inact."),
                  paste(gene1,"act.,",gene2,"act.")
                  ),col=c("darkblue","forestgreen","black","darkred"),pch=c(16,15,4,17),cex=font_cex,bty="n")
  #box("plot", col="red")
  #box("figure", col="blue")
  
  mtext(paste("Expression distribution for",gene1,"and",gene2), adj=0.5, side=3, outer=TRUE, cex=font_cex)
}

### Format:
### compare_genes(transcription_factor, target, transcription_factor_name, target_name)
#pdf(file="/home/jason/Documents/School/College/Fall Junior/STAT RESEARCH/crp", width=13, height=8)
# compare_genes(65,63,'araC','araA',3292,"crp")
# compare_genes(65,64,'araC','araB',3292,"crp")
# compare_genes(65,62,'araC','araD',3292,"crp")
#dev.off()
#compare_genes(2643,459,'mprA','acrA')
#compare_genes(3956,293,'zur','rpmE2')
#compare_genes(1323,725,'fnr','cydA')
#compare_genes(502,501,'allR','allA')

# #pdf(file="/home/jason/Documents/School/College/Fall Junior/STAT RESEARCH/TRN Conflicts UniMM ", width=13, height=8)
# compare_genes(name_map,'araC','araA','crp')
# compare_genes(name_map,'araC','araB','crp')
# compare_genes(name_map,'araC','araD','crp')
# compare_genes(name_map,'araC','araF','crp')
# compare_genes(name_map,'araC','araG','crp')
# #compare_genes(name_map,'araC','araH','crp') no peg id
# compare_genes(name_map,'araC','ytfQ','crp')
# #compare_genes(name_map,'araC','ytfR','crp') no peg id
# compare_genes(name_map,'araC','ytfT','crp')
# compare_genes(name_map,'araC','yjfF','crp') 
# #dev.off()

plot.elaborate("araC","araD",0.5,0.5)
