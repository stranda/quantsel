###summarise a simulation replicated in a list
###
simsum <- function(allreslst, CORES=1, plt=F)
{
    treatdf <- as.data.frame(do.call(rbind,lapply(allreslst,function(x){x[[1]][3][[1]]})))
    allreslst <- allreslst[!is.na(treatdf[,1])]
    treatdf <- treatdf[!is.na(treatdf[,1]),]
    treatdf$rep <- 1:nrow(treatdf)    

###calculate the size of the dataframe required
    hn <- nrow(allreslst[[1]][[1]]$phenosum)
    gn <- length(allreslst[[1]])
    pn <- length(grep("mean",names(allreslst[[1]][[1]]$phenosum)))
    rp <- length(allreslst)
    rows <- hn*gn*pn*rp
    allPhensum <- as.data.frame(matrix(NA,nrow=rows,ncol=9))
    names(allPhensum) <- c("gen","pop","n","rep","phenotype","mean","seedmix","repro","denstol")    
    for (r in 1:length(allreslst))
    {
        phenosumdf <- do.call(rbind,lapply(allreslst[[r]],function(x){x$phenosum}))
        phenosumdf <- phenosumdf[,!grepl("_sd",names(phenosumdf))]
        phenosumdf$rep <- r
        phenosumdf <- pivot_longer(phenosumdf,cols=ends_with(c("mean"))) %>%
            mutate(phenotype=gsub("_mean","",gsub("phen","",name)),mean=value) %>%
            select(-name,-value)
        phenosumdf <- phenosumdf %>% left_join(treatdf)
        allPhensum[((r-1)*rows+1):((r-1)*rows+rows),] <- as.data.frame(phenosumdf)
    }
    
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
