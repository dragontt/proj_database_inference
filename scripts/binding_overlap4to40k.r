args <- commandArgs(trailingOnly = TRUE)
network <- toString(args[1])
network.name <- toString(args[2])

net <- read.table(network,header=TRUE)
bsinfo <- read.table("~zmaier/Yeast/PWM/PBMs/binding_potential/binding_summary/all.interactions.summary.combined.orfs.named",header=TRUE)
pdna.inter <- read.table("~bhaynes/projects/network_inference/figures/yeast_meta/pwm_analysis/binding_site_analysis_chip/chip.interactions",header=TRUE)
source("~bhaynes/projects/network_inference/util/util.r")

chip.bp.np.setsizes <- c()
chip.bp.strong.np.setsizes <- c()
chip.bp.weak.np.setsizes <- c()


chip_desired_recall <- 0.33
np_conf_cutoff <- 0.02036

orf_universe <- intersect(union(unique(bsinfo$TARGET),unique(pdna.inter$TARGET)),unique(net$TARGET))
reg_universe <- intersect(unique(as.character(pdna.inter$REGULATOR)),unique(as.character(bsinfo$REGULATOR)))
interaction_universe_cnt <- length(reg_universe) * length(orf_universe)

np_conf_cutoffs <- net$CONFIDENCE[order(net$CONFIDENCE,decreasing=TRUE)][seq(4000,40000,by=4000)]
for (np_conf_cutoff in np_conf_cutoffs)
{
	for (chip_desired_recall in c(0.1,0.33,0.5))
	{
		bp.cutoff <- c()
		pdna.target.count <- c()
		predicted.bp.support <- c()
		
		chip.bp.np.sets <- c()
		chip.bp.strong.np.sets <- c()
		chip.bp.weak.np.sets <- c()

		i <- 1
		for (regulator in intersect(unique(as.character(pdna.inter$REGULATOR)),unique(as.character(bsinfo$REGULATOR))) )
		{
			targets <- union(as.character(pdna.inter$TARGET[which(pdna.inter$REGULATOR==regulator)]),as.character(bsinfo$TARGET[which(bsinfo$REGULATOR==regulator)]))
			pdna.evid <- rep(0,times=length(targets))
			binding.evid.max <- rep(0,times=length(targets))
			binding.evid.sum <- rep(0,times=length(targets))
			pdna.evid[match(as.character(pdna.inter$TARGET[which(pdna.inter$REGULATOR==regulator)]),targets)] <- 1
			binding.evid.max[match(as.character(bsinfo$TARGET[which(bsinfo$REGULATOR==regulator)]),targets)] <- bsinfo$MAXP[which(bsinfo$REGULATOR==regulator)]
			binding.evid.sum[match(as.character(bsinfo$TARGET[which(bsinfo$REGULATOR==regulator)]),targets)] <- bsinfo$SUMP[which(bsinfo$REGULATOR==regulator)]
			cat(regulator,length(targets),"\n")
			pdna.target.count[i] <- sum(pdna.evid)
			binding.evid.comb <- apply(rbind(binding.evid.max,binding.evid.sum),2,max)	
			prc.max <- compute.prc(binding.evid.max,pdna.evid)
			prc.sum <- compute.prc(binding.evid.sum,pdna.evid)
			prc.comb <- compute.prc(binding.evid.comb,pdna.evid)
		
			if 	(length(which(prc.comb[,1]<=chip_desired_recall))>0)
			{
				bp.cutoff[i] <- binding.evid.comb[sort.list(binding.evid.comb,decreasing=TRUE)][length(which(prc.comb[,1]<=chip_desired_recall))]
			}
			else
			{
				bp.cutoff[i] <- 1
			}
		
			predicted.targets <- as.character(net$TARGET[which(net$REGULATOR==regulator)])[which(net$CONFIDENCE[which(net$REGULATOR==regulator)]>np_conf_cutoff)]
			predicted.targets <- intersect(predicted.targets,orf_universe) # Narrow to space under consideration

			predicted.targets.chip.supported <- intersect(predicted.targets,as.character(pdna.inter$TARGET[which(pdna.inter$REGULATOR==regulator)]))
			predicted.targets.chip.unsupported <- setdiff(predicted.targets,as.character(pdna.inter$TARGET[which(pdna.inter$REGULATOR==regulator)]))

			bp.targets <- targets[which(binding.evid.comb>bp.cutoff[i])]
			bp.strong.targets <- targets[which(binding.evid.max>bp.cutoff[i])]
			bp.weak.targets <- setdiff(targets[which(binding.evid.sum>bp.cutoff[i])],bp.strong.targets)

			pdna.targets <- targets[which(pdna.evid>0)]
			
			chip.bp.np.sets <- rbind(chip.bp.np.sets,c(length(bp.targets),sum(pdna.evid),length(predicted.targets),length(intersect(bp.targets,pdna.targets)),length(intersect(bp.targets,predicted.targets)),length(intersect(predicted.targets,pdna.targets)),length(intersect(bp.targets,intersect(predicted.targets,pdna.targets)))))

			chip.bp.strong.np.sets <- rbind(chip.bp.strong.np.sets,c(length(bp.strong.targets),sum(pdna.evid),length(predicted.targets),length(intersect(bp.strong.targets,pdna.targets)),length(intersect(bp.strong.targets,predicted.targets)),length(intersect(predicted.targets,pdna.targets)),length(intersect(bp.strong.targets,intersect(predicted.targets,pdna.targets)))))

				chip.bp.weak.np.sets <- rbind(chip.bp.weak.np.sets,c(length(bp.weak.targets),sum(pdna.evid),length(predicted.targets),length(intersect(bp.weak.targets,pdna.targets)),length(intersect(bp.weak.targets,predicted.targets)),length(intersect(predicted.targets,pdna.targets)),length(intersect(bp.weak.targets,intersect(predicted.targets,pdna.targets)))))
	
			i <- i + 1
		}
		
		chip.bp.np.setsize <- apply(chip.bp.np.sets,2,sum)
		chip.bp.strong.np.setsize <- apply(chip.bp.strong.np.sets,2,sum)
		chip.bp.weak.np.setsize <- apply(chip.bp.weak.np.sets,2,sum)
		
		chip.bp.np.setsizes <- rbind(chip.bp.np.setsizes,chip.bp.np.setsize)
		chip.bp.strong.np.setsizes <- rbind(chip.bp.strong.np.setsizes,chip.bp.strong.np.setsize)
		chip.bp.weak.np.setsizes <- rbind(chip.bp.weak.np.setsizes,chip.bp.weak.np.setsize)
	
	}
}

write.table(chip.bp.np.setsizes,file="chip.bp.np.set.sizes.top4to40k.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')

## Make NetPropet Figures
results <- chip.bp.np.setsizes
pwm_random <- results[1,1]/interaction_universe_cnt
chip_random <- results[1,2]/interaction_universe_cnt
cat(pwm_random,"\t", chip_random,"\n")
num_confidence_intervals <- dim(results)[1]/3
pwm_overlaps <- vector(mode="numeric", length=num_confidence_intervals)
chip_overlaps <- vector(mode="numeric", length=num_confidence_intervals)
for(i in 1:num_confidence_intervals) {
	lower_confidence <- (i-1) * 4000
	upper_confidence <- i * 4000
	row <- (3* i) -2
	if(i ==1 ) {
		pwm_overlaps[i] <-  results[row,5]/results[row,3]
		chip_overlaps[i] <-  results[row,6]/results[row,3]
		strong_combined <- c(results[row,5]/results[row,3], results[row,4]/results[row,2], results[row,7]/results[row,6])
		row <- row + 1
		medium_combined <- c(results[row,5]/results[row,3], results[row,4]/results[row,2], results[row,7]/results[row,6])
		row <- row + 1
		weak_combined <- c(results[row,5]/results[row,3], results[row,4]/results[row,2], results[row,7]/results[row,6])
		
	} else {
		pwm_overlaps[i] <-  (results[row,5] - results[row-3,5])/(results[row,3] - results[row-3,3])
		chip_overlaps[i] <-  (results[row,6] - results[row-3,6])/(results[row,3] - results[row-3,3])
		strong_combined <- c((results[row,5] - results[row-3,5])/(results[row,3] - results[row-3,3]), results[row,4]/results[row,2], (results[row,7] - results[row-3,7])/(results[row,6] - results[row-3,6]))
		row <- row + 1
		medium_combined <- c((results[row,5] - results[row-3,5])/(results[row,3] - results[row-3,3]), results[row,4]/results[row,2], (results[row,7] - results[row-3,7])/(results[row,6] - results[row-3,6]))
		row <- row + 1
		weak_combined <- c((results[row,5] - results[row-3,5])/(results[row,3] - results[row-3,3]), results[row,4]/results[row,2], (results[row,7] - results[row-3,7])/(results[row,6] - results[row-3,6]))
	}
	if(i ==1) {
                outname = paste("ChIP_Support_PWM_vs_NP_BarPlot",
                                network.name, ".jpg", sep="")
		jpeg(outname)
		barplot(cbind(strong_combined, medium_combined, weak_combined), beside=TRUE, col=c(2,4,3), names.arg=c("High", "Medium", "Low"), legend.text=c("NetProphet", "ChIP", "ChIP & NetProphet"), main=paste("Binding Potential Analysis\nNetProphet Scores: ", lower_confidence, " - ", upper_confidence, sep=""), xlab="Stringency of threshold for PWM support", ylab="Ineractions supported by PWM models", ylim=c(0,0.8))
		dev.off()
	}
}

outname = paste("PWM_Support_On_NP_Confidence_Bins_Line_Plot",
                network.name, ".jpg", sep="")
jpeg(outname)
plot(seq(4000,4000*num_confidence_intervals, 4000), pwm_overlaps, col=3, xlab="Predictions grouped by rank", ylab="Interactions supported by PWMs", xaxp =c(4000, 4000*num_confidence_intervals, num_confidence_intervals -1 ), ylim=c(0,1.2*max(pwm_overlaps)), yaxp=c(0,0.18,6))
lines(seq(4000,4000*num_confidence_intervals, 4000), pwm_overlaps, col=3)
lines(seq(4000,4000*num_confidence_intervals, 4000), rep(pwm_random,num_confidence_intervals), lty=2, col="light green")
dev.off()

outname = paste("ChIP_Support_On_NP_Confidence_Bins_Line_Plot",
                network.name, ".jpg", sep="")
jpeg(outname)
plot(seq(4000,4000*num_confidence_intervals, 4000), chip_overlaps, col=4, xlab="Predictions grouped by rank", ylab="Interactions supported by ChIP", xaxp =c(4000, 4000*num_confidence_intervals, num_confidence_intervals -1 ), ylim=c(0,1.2*max(chip_overlaps)), yaxp=c(0,0.25,5))
lines(seq(4000,4000*num_confidence_intervals, 4000), chip_overlaps, col=4)
lines(seq(4000,4000*num_confidence_intervals, 4000), rep(chip_random,num_confidence_intervals), lty=2, col="light blue")
dev.off()

pwm_overlaps_new <- pwm_overlaps
chip_overlaps_new <- chip_overlaps

results <- as.matrix(read.table("../netprophet_trained/chip.bp.np.set.sizes.top4to40k.txt"))
pwm_overlaps <- vector(mode="numeric", length=num_confidence_intervals)
chip_overlaps <- vector(mode="numeric", length=num_confidence_intervals)
for(i in 1:num_confidence_intervals) {
	lower_confidence <- (i-1) * 4000
	upper_confidence <- i * 4000
	row <- (3* i) -2
	if(i ==1 ) {
		pwm_overlaps[i] <-  results[row,5]/results[row,3]
		chip_overlaps[i] <-  results[row,6]/results[row,3]
		strong_combined <- c(results[row,5]/results[row,3], results[row,4]/results[row,2], results[row,7]/results[row,6])
		row <- row + 1
		medium_combined <- c(results[row,5]/results[row,3], results[row,4]/results[row,2], results[row,7]/results[row,6])
		row <- row + 1
		weak_combined <- c(results[row,5]/results[row,3], results[row,4]/results[row,2], results[row,7]/results[row,6])
		
	} else {
		pwm_overlaps[i] <-  (results[row,5] - results[row-3,5])/(results[row,3] - results[row-3,3])
		chip_overlaps[i] <-  (results[row,6] - results[row-3,6])/(results[row,3] - results[row-3,3])
		strong_combined <- c((results[row,5] - results[row-3,5])/(results[row,3] - results[row-3,3]), results[row,4]/results[row,2], (results[row,7] - results[row-3,7])/(results[row,6] - results[row-3,6]))
		row <- row + 1
		medium_combined <- c((results[row,5] - results[row-3,5])/(results[row,3] - results[row-3,3]), results[row,4]/results[row,2], (results[row,7] - results[row-3,7])/(results[row,6] - results[row-3,6]))
		row <- row + 1
		weak_combined <- c((results[row,5] - results[row-3,5])/(results[row,3] - results[row-3,3]), results[row,4]/results[row,2], (results[row,7] - results[row-3,7])/(results[row,6] - results[row-3,6]))
	}
}

outname = paste("PWM_Support_On_NP_Confidence_Bins_Line_Plot_Comparison_",
                network.name, ".jpg", sep="")
jpeg(outname)
plot(NULL, col=3, xlab="Predictions grouped by rank", ylab="Interactions supported by PWMs", xlim=c(4000,4000*num_confidence_intervals), xaxp =c(4000, 4000*num_confidence_intervals, num_confidence_intervals -1 ), ylim=c(0,1.2*max(pmax(pwm_overlaps, pwm_overlaps_new))), yaxp=c(0,0.18,6))
lines(seq(4000,4000*num_confidence_intervals, 4000), pwm_overlaps, col=3)
lines(seq(4000,4000*num_confidence_intervals, 4000), pwm_overlaps_new, col=1)
lines(seq(4000,4000*num_confidence_intervals, 4000), rep(pwm_random,num_confidence_intervals), lty=2, col="light green")
legend(0.6*(4000*num_confidence_intervals), 0.16, legend=c("NetProphet", "Improvement?"), lwd=1, col=c(3,1))
dev.off()

outname = paste("ChIP_Support_On_NP_Confidence_Bins_Line_Plot_Comparison",
                network.name, ".jpg", sep="")
jpeg(outname)
plot(NULL, col=4, xlab="Predictions grouped by rank", ylab="Interactions supported by ChIP", xlim=c(4000,4000*num_confidence_intervals), xaxp =c(4000, 4000*num_confidence_intervals, num_confidence_intervals -1 ), ylim=c(0,1.2*max(pmax(chip_overlaps, chip_overlaps_new))), yaxp=c(0,0.25,5))
lines(seq(4000,4000*num_confidence_intervals, 4000), chip_overlaps, col=4)
lines(seq(4000,4000*num_confidence_intervals, 4000), chip_overlaps_new, col=1)
lines(seq(4000,4000*num_confidence_intervals, 4000), rep(chip_random,num_confidence_intervals), lty=2, col="light blue")
legend(0.6*(4000*num_confidence_intervals), 0.25, legend=c("NetProphet", "Improvement?"), lwd=1, col=c(4,1))
dev.off()

