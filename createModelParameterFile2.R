createModelParameterFile2 <- function(startyear=2000,endyear=as.numeric(substr(Sys.Date(),1,4)),hmin=80,hmax=150,hres=2,fname='ElSpecModelParameters'){
    #
    #
    # Tabulate the necessary IRI parameters in Tromso and Svalbard
    # We will need O2+, NO+, and Te, but also O+ and Ti are tabulated
    # because they may be very useful outside ELSPEC
    #
    #
    #
    #
    # IV 2017, 2022
    #
    # Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
    # This is free software, licensed under GNU GPL version 2 or later


    require(IRI2016)
    require(R.matlab)

    years <- seq(startyear,endyear)
    nyear <- length(years)

    h <- seq(hmin,hmax,by=hres)
    nh <- length(h)

    ntot <- nyear*12*31*24

    ns <- 0
    for(iyear in seq(nyear)){
        iriparTRO<-array(dim=c(12,31,24,nh,5))
        iriparESR<-array(dim=c(12,31,24,nh,5))
        for(month in seq(1,12)){
            for(day in seq(1,31)){
                for(hour in seq(0,23)){
                    ns <- ns+1
                    if (years[iyear]>1980){
                        tmpTRO<-tryCatch(iriParams(time=c(years[iyear],month,day,hour,30,0),heights=h,latitude=69.5864,longitude=19.2272), error=function(e){return(NULL)});
                        if (!is.null(tmpTRO)){
                            iriparTRO[month,day,hour+1,,]<-t(tmpTRO[c('O+','O2+','NO+','Ti','Te'),])
                        }
                    }
                    if(years[iyear]>1995){
                        tmpESR<-tryCatch(iriParams(time=c(years[iyear],month,day,hour,30,0),heights=h,latitude=78.15,longitude=16.02), error=function(e){return(NULL)});
                        if (!is.null(tmpESR)){
                            iriparESR[month,day,hour+1,,]<-t(tmpESR[c('O+','O2+','NO+','Ti','Te'),])
                        }
                    }
                    cat('                      \r',sprintf("%6.0f %2.0f %2.0f %2.0f %3.0f",years[iyear],month,day,hour,ns/ntot*100,'%'))
                }
            }
        }
        if(years[iyear]>1980){
            writeMat(paste(fname,'TRO',as.character(years[iyear]),'.mat',sep=''),iripar=iriparTRO)
        }
        if(years[iyear]>1995){
            writeMat(paste(fname,'ESR',as.character(years[iyear]),'.mat',sep=''),iripar=iriparESR)
        }
    }

}
