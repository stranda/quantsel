###summarise a simulation replicated in a list
###
simsum <- function(allreslst, fn=paste0("test_res.rda"))
{
    

allPhensum <- do.call(rbind,mclapply(1:length(allreslst),mc.cores=CORES,function(ar)
{
    phenosumdf <- do.call(rbind,lapply(allreslst[[ar]],function(x){x$phenosum}))
    phenosumdf <- phenosumdf[,!grepl("_sd",names(phenosumdf))]
    phenosumdf$rep=ar
    phenosumdf <- pivot_longer(phenosumdf,cols=ends_with(c("mean"))) %>%
        mutate(phenotype=gsub("_mean","",gsub("phen","",name)),mean=value) %>% select(-name,-value)
    phenosumdf
}))

allAfrqsum <- do.call(rbind,mclapply(1:length(allreslst),mc.cores=CORES,function(ar)
{
    afreqsumdf <- do.call(rbind,lapply(allreslst[[ar]],function(x){x$afreq}))
    afreqsumdf$rep=ar
    afreqsumdf
}))

save(file=fn,allreslst,allPhensum,allAfrqsum)

print(
allPhensum %>%
    ggplot(aes(x=gen,y=mean,color=phenotype)) +
    geom_smooth() +
    facet_wrap (~pop)
)
    
}
