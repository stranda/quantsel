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
        gens <- c(allgens)
    } else if (gen==0)
    {
        gens <- range(allgens)[1]
    } else gens <- range(allgens)[2]
    
    focal_db <- phen_db %>%
        filter(gen %in% gens) %>%
        select(all_of(onm),all_of(phennm))
    if (length(gens)==1)
        {
            for (p in phennm)
            {
                ggplot(focal_db,aes(x=pop,y=get(p),group=seedmix,color=seedmix)) +geom_smooth() +
                    facet_wrap(~repro+denstol,labeller=label_both) + ggtitle(paste("Phenotype: ",p," Generation: ",gens))+
                    ylim(c(0,1))+
                    ylab(p) + xlab("transect coordinate 0=creek")
            }
        }
}
