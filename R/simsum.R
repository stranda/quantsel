###summarise a simulation replicated in a list
###
simsum <- function(allreslst, CORES=1, plt=F)
{
    treatdf <- as.data.frame(do.call(rbind,lapply(allreslst,function(x){cbind(x[[1]]$treat)})))
    treatdf$rep <- 1:nrow(treatdf)    

    allPhensum <- do.call(rbind,mclapply(1:length(allreslst),mc.cores=CORES,function(ar)
    {
        phenosumdf <- do.call(rbind,lapply(allreslst[[ar]],function(x){x$phenosum}))
        phenosumdf <- phenosumdf[,!grepl("_sd",names(phenosumdf))]
        phenosumdf$rep=ar
        phenosumdf <- pivot_longer(phenosumdf,cols=ends_with(c("mean"))) %>%
            mutate(phenotype=gsub("_mean","",gsub("phen","",name)),mean=value) %>%
            select(-name,-value)
        phenosumdf %>% left_join(treatdf)
    }))
    
    allAfrqsum <- do.call(rbind,mclapply(1:length(allreslst),mc.cores=CORES,function(ar)
    {
        afreqsumdf <- do.call(rbind,lapply(allreslst[[ar]],function(x){x$afreq}))
        afreqsumdf$rep=ar
        afreqsumdf %>% left_join(treatdf)
    }))

    
    if (plt==T)
        print(
            allPhensum %>%
            ggplot(aes(x=gen,y=mean,color=phenotype)) +
            geom_smooth() +
            facet_wrap (~pop)
        )
    list(phen=allPhensum,af=allAfrqsum)
}
