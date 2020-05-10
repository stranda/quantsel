###
### summarise phenotypes across one dimension 
###
###
#' @export
plotPhenos1d <- function(dbfile="results.db",pheno=NULL,
                         gen=1 #if gen=0 plot first gen, if 1 plot last gen, if NULL, try to plot all
                        )

{
    if (!file.exists(dbfile)) stop (paste("Can't find ",dbfile))
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = dbfile) #connect to sqlite db on disk
    phen_db <- tbl(con,"phen")
    nm <- names(RSQLite::dbGetQuery(con,"SELECT * FROM phen LIMIT 1"))
    phennm <- nm[grep("phen",nm)]
    onm <- nm[-grep("phen",nm)]
    if (!is.null(pheno))
    {
        nums <- as.numeric(gsub("_mean|_sd","",gsub("phen","",phennm)))
        phennm <- phennm[which(nums %in% pheno)]
    }
    phennm <- phennm[grep("mean",phennm)]
    allgens <-  collect(select(phen_db,gen) %>% distinct() )
    if (is.null(gen))
    {
        focal_db <- phen_db %>%
            select(all_of(onm),all_of(phennm))
    } else if (gen==0)
    {
        gens <- range(allgens)[1]
        focal_db <- phen_db %>% filter(gen %in% gens) %>%
            select(all_of(onm),all_of(phennm))
    } else {
        gens <- range(allgens)[2]
        focal_db <- phen_db %>% filter(gen %in% gens) %>%
            select(all_of(onm),all_of(phennm))
    }
    focal_db <- mutate(focal_db,zone=ifelse(pop<8,"tall",ifelse(pop>16,"short","intermediate")))
#    plotlst <- vector("list",length(phennm))
#    names(plotlst) <- sapply(phennm,function(p) {paste0(p,"_trans")})
    
    if (!is.null(gen))
    {
            for (p in phennm)
            {
                plt <- ggplot(focal_db,aes(x=pop,y=get(p),group=seedmix,color=seedmix)) + geom_smooth() +
                    facet_wrap(~repro+denstol,labeller=label_both) + ggtitle(paste("Phenotype: ",p," Generation: ",gens))+
                    ylim(c(0,1))+
                    ylab(p) + xlab("transect coordinate 0=creek")

                print(plt)
            }
    } else {
        for (p in phennm)
        {
            plt <- ggplot(focal_db,aes(x=gen,y=get(p),group=seedmix,color=seedmix)) + geom_smooth() +
                facet_wrap(~repro+denstol,labeller=label_both) + ggtitle(paste("Phenotype: ",p))+
                ylim(c(0,1))+
                ylab(p) + xlab("Time click")
            print(plt)

            plt <- ggplot(focal_db,aes(x=gen,y=get(p),group=seedmix,color=seedmix)) + geom_smooth() +
                facet_wrap(repro+denstol~zone,labeller=label_both) + ggtitle(paste("Phenotype: ",p))+
                ylim(c(0,1))+
                ylab(p) + xlab("Time click")
            print(plt)

            for (sm in unique(collect(focal_db)$seedmix))
                {
                    surf_db <- filter(focal_db,seedmix==sm) %>% group_by(pop,gen,denstol,repro) %>%
                        summarise_all(.funs=mean)
                    plt <- ggplot(surf_db,aes(y=gen,x=pop,z=get(p))) + geom_contour(bin=4) +
                        facet_wrap(repro~denstol,labeller=label_both)+ggtitle(paste("Phenotype: ",p,"Seedmix: ",sm))
                    
                    print(plt)
                }
        }
        
    }
}
