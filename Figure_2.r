coverage_sense=read.csv("~/Downloads/OneDrive_2025-07-16/Analyse_Camille/HVT/aln_sense_coverage.csv",sep="\t",header=FALSE)
coverage_antisense=read.csv("~/Downloads/OneDrive_2025-07-16/Analyse_Camille/HVT/aln_antisense_coverage.csv",sep="\t",header=FALSE)
coverage_virus=read.csv("~/Downloads/OneDrive_2025-07-16/Analyse_Camille/HVT/aln_virus_coverage.csv",sep="\t",header=FALSE)
head(coverage_virus)
totcov=sum(coverage_virus$V3)
head(coverage_sense)
ggplot(coverage_sense,aes(x=V2,y=(V3/totcov)*100000000)) + geom_line() + geom_line(data=coverage_antisense,aes(x=V2,y=(-V3/totcov)*100000000)) + theme_minimal() + xlab("Position on the genome") + ylab("Coverage (RPB)") + ggtitle("Coverage of the sense and antisense reads") + xlim(123435,136989)+ylim(-30,250)
ggsave("~/Downloads/OneDrive_2025-07-16/Analyse_Camille/HVT/coverage_sense_antisense.png",width=12,height=4)






