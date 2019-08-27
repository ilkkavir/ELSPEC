createModelParameterFile <- function(startyear=2000,endyear=as.numeric(substr(Sys.Date(),1,4)),tres=1,hmin=80,hmax=150,hres=2,fname='ElSpecModelParameters'){
    #
    #
    # Tabulate necesary IRI and MSIS parameters and save them in a matlab mat-file
    #
    # IV 2017
    #
    # Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
    # This is free software, licensed under GNU GPL version 2 or later


    require(IRI2016)
    require(R.matlab)

    years <- seq(startyear,endyear)
    nyear <- length(years)

    h <- seq(hmin,hmax,by=hres)
    nh <- length(h)

    t <- seq(1,24,by=tres)
    nt <- length(t)


    ntot <- nyear*12*31*nt

    ns <- 0
    for(iyear in seq(nyear)){
        iripar<-array(dim=c(12,31,nt,nh,10))
        for(month in seq(1,12)){
            for(day in seq(1,31)){
                for(hour in seq(1,24)){
                    ns <- ns+1
                    tmp<-tryCatch(iriParams(time=c(years[iyear],month,day,hour,0,0),heights=h),
                                  error=function(e){return(NULL)});
                    if (!is.null(tmp)){
                        iripar[month,day,hour,,]<-t(tmp[c('O+','O2+','NO+','cluster','O','N2','O2','Tn','Ti','Te'),])
                    }
                    cat('                      \r',sprintf("%6.0f %2.0f %2.0f %2.0f %3.0f",years[iyear],month,day,hour,ns/ntot*100,'%'))
                }
            }
        }
        writeMat(paste(fname,as.character(years[iyear]),'.mat',sep=''),iripar=iripar)
    }

}
