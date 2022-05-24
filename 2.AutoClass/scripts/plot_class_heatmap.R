## (c) Mikhail Spivakov 2011

library(gplots)


cbReadPlotSave = function(class_tab_fname, plot_fname, qbm, ct_idname="V1", qbm_idname="CRMID", merge=vector("list"), classcol=3, ...){
	
	cqbm = make_cqbm(class_tab_fname, qbm, ct_idname, qbm_idname)

	if (length(merge)){
		
		for (i in 1:length(merge)){
			start = merge[[i]][1]
			for (j in 2:length(merge[[i]])){
				cat ("merging class ", merge[[i]][j], " with ", merge[[i]][1], "\n")
				cqbm[cqbm[,classcol]==merge[[i]][j],classcol] = merge[[i]][1]
			}
		}

	}
	

	plot_class_barcharts(cqbm, classcol, ...) 
	readline("Choose a nice window shape and press any key...")
	quartz.save(plot_fname, type="pdf")
	
	invisible(cqbm)

} 


cbReadPlotClassSave = function(class_tab_fname, whichclass, qbm, qbm_idcol=1, classcol=3, crmidcol=1, qbm_idname="CRMID", ct_idname="V1", plot_idname="%qbmid%", plot_fname,  ...){
	a = make_cqbm(class_tab_fname, qbm, ct_idname, qbm_idname)
	cbPlotEachSave(a[a[,classcol]==whichclass, crmidcol],plot_fname, qbm=qbm, qbm_idcol=qbm_idcol, plot_idname=plot_idname, ...)
}

cbPlotEachSave = function(crmlist, plot_fname, qbm, qbm_idcol=1, tfstcol = 5, tfendcol = 14, min=-1, max=-1, plot_idname="%qbmid%",  ...){
	
	if (length(crmlist)==0){
		cat("Error: CRM list is empty\n")
		return(-1)	
	}
	
	cqbm = qbm[ qbm[,qbm_idcol] %in% crmlist , ]
	
	if (plot_idname=="%qbmid%"){
		plot_idname = names(cqbm)[qbm_idcol]
	}			
	
	if ( min > 0){
		cqbm = cqbm[ min:length(cqbm[,1]), ]
	}
	
	if (max > 0){
		cqbm = cqbm[ 1:max, ]	
	}
		
	plot_class_barcharts(cqbm, classcol = which(names(cqbm)==plot_idname), tfstcol=tfstcol, tfendcol=tfendcol, single=T, ...)
	readline("Choose a nice window shape and press any key...")
	quartz.save(plot_fname, type="pdf")
	
		
}

chReadPlotSave = function(class_tab_fname, plot_fname, qbm, ct_idname="V1", qbm_idname="CRMID", merge=vector("list"), classcol=3, ...){

	cqbm = make_cqbm(class_tab_fname, qbm, ct_idname, qbm_idname)

	if (length(merge)){
		
		for (i in 1:length(merge)){
			start = merge[[i]][1]
			for (j in 2:length(merge[[i]])){
				cat ("merging class ", merge[[i]][j], " with ", merge[[i]][1], "\n")
				cqbm[cqbm[,classcol]==merge[[i]][j],classcol] = merge[[i]][1]
			}
		}

	}
	

	plot_class_heatmap(cqbm, classcol, ...) 
	readline("Choose a nice window shape and press any key...")
	quartz.save(plot_fname, type="pdf")
}

make_cqbm = function (class_tab_fname, qbm, ct_idname="V1", qbm_idname="CRMID"){
	a = read.table(class_tab_fname)
	cqbm = merge(a,qbm,by.x=ct_idname,by.y=qbm_idname)
	cqbm
}


plot_class_heatmap = function(cqbm,classcol=3,tfstcol=7,tfendcol=16,rel=T, ntimes=NULL, pal=greenred(100), verbose=F, 
                              scalebar=F, sym=T, zlim=T, per.row=4, single=F){

	#zlim is deprecated - sym is a better name for what it's now doing
	if (zlim == F){
		sym = F
		cat ("Warning: zlim is deprecated. Use sym instead next time\n")
	}
	
	cat("\n")
	
	res = vector("list")
	rms = vector("list")
	
	nms = names(cqbm)[tfstcol:tfendcol]

	tfs = gsub("(.+)_.+","\\1",nms)
	tfs = unique(tfs)
	ntfs = length(tfs)

	if (!is.factor(cqbm[,classcol])){
		cqbm[,classcol] = as.factor(cqbm[,classcol])
	}
	
	else{
		if (is.integer(levels(cqbm[,classcol]))){
			cqbm[,classcol] = as.factor(as.integer(cqbm[,classcol]))
		}
		else{
			cqbm[,classcol] = as.factor(as.character(cqbm[,classcol]))
		}
	}
	#this refactoring is done because there may be unused levels 
	#in the input cqbm, as shown by experience
		
	if (verbose==T){
		print(tfs)
	}
	
  if(is.null(ntimes)){
  	times = gsub(".+_(.+)","\\1",nms) 
	  times = unique(times)
	  ntimes = length(times)
  }
  else{
    times=""
  }
  

	if (verbose==T){
		cat("Calculating total means...\n")
	}
	
	tm = matrix(data=numeric(0),nrow=ntimes,ncol=ntfs)
	for (j in 1:ntfs){
		for (i in 1:ntimes){
			tm[i,j] = mean(cqbm[,tfstcol+(j-1)*ntimes+i-1])
#			cat("\t",names(cqbm)[tfstcol+(j-1)*ntimes+i-1],"\n")
		}
	}
	
	if (verbose==T){
		print(tm)
		cat("\n")
	}

#	maxc = length(levels(cqbm[,classcol]))
#	print(paste("maxc=",maxc))
	
	sb=0
	if (scalebar==T){
		sb=1
	}
	
	rows = ceiling((length(levels(cqbm[,classcol]))+sb)/per.row)
	par(mfrow=c(rows,per.row))
	
	tmin = 10000000
	tmax = -10000000
				
	for(cc in levels(cqbm[,classcol])){
		
		
		if (verbose==T){ 
			cat("\nProcessing class",cc,"...\n")
		}
		
		m = matrix(data=numeric(0),nrow=ntimes,ncol=ntfs)
		rm = matrix(data=numeric(0),nrow=ntimes,ncol=ntfs)

		for (j in 1:ntfs){
			for (i in 1:ntimes){
				m[i,j] = mean(cqbm[cqbm[,classcol]==cc,tfstcol+(j-1)*ntimes+i-1])
#				cat("\t",names(cqbm)[tfstcol+(j-1)*ntimes+i-1],"\n")
			}	
		}


		if (rel == T){
				rm = m-tm
		}
		else{
				rm = m
		}
		
		
		res[[as.character(cc)]] = m	
		rms[[as.character(cc)]] = rm

		if (verbose==T){
			cat ("min rm=", min(rm), "tmin=", tmin, "\n")
		}
		
		if (min(rm)<tmin){
			tmin <- min(rm)
		}
		if (max(rm)>tmax){
			tmax <- max(rm)
		}
		
	}


	if (sym==T){
		if (abs(tmin)>abs(tmax)){
			zlimc = c(tmin,-1*tmin)
		}
		else{
			zlimc = c(-1*tmax, tmax)
		}
	}
	else{
		zlimc = c(tmin, tmax)
	}
	
	if (verbose == T){
		print (rms)
		cat ("total min = ", tmin, "total max = ", tmax, "zlimc = {", zlimc[1], ", ", zlimc[2], "}\n")
	}

	for (cc in levels (cqbm[,classcol])){
		
		x = 100*(1:ncol(rms[[as.character(cc)]]))
		y = 200*(1:nrow(rms[[as.character(cc)]]))

		image(x,y,t(rms[[as.character(cc)]]),zlim=zlimc, col=pal,axes=F, xlab="", ylab= "")
		if (single==T){
			title (cc)
		}
		else{
			title (paste("class",cc,"( n =",length(cqbm[cqbm[,classcol]==cc,1]),")"))
		}
		axis(1,at=x,labels=tfs)	
		axis(2,at=y,labels=times)		
	}
	
	if (sb==1){
		sb = matrix(seq(min(zlimc), max(zlimc), length.out=20),c(1,20))
		y = 100*(1:nrow(sb))
		x = 100*(1:ncol(sb))
		image(x,y,t(sb),zlim=zlimc, col=pal, axes=F,xlab="", ylab="") 
		axis(1,at=x,labels=apply(sb,2, function(x)sprintf("%.2f",x)))
		title ("scalebar")	
	}

	res[["total"]] = tm

	if (verbose==T){
		print (paste("range res:", range (res)))
	}

	
	invisible(res)
}

#plot_class_heatmap(aq05)

plot_class_barcharts = function(cqbm,classcol=3,tfstcol=7,tfendcol=16, tfrange=NULL, 
                                upperonly=T, legend.text=F, plot.ci=T, per.row=4, single=F, 
                                commonylim=TRUE, refactor=TRUE, ylim=NULL, ...){
	
	res = vector("list")
		
	nms = NULL
	contrange = TRUE

	if (is.null(tfrange)){
	  tfrange = tfstcol:tfendcol
	}
	else{
		tfstcol = min(tfrange)
		tfendcol = max(tfrange)
	}

	nms = names(cqbm)[tfrange]	
	  
	if(! all(min(tfrange):max(tfrange) %in% tfrange)){
			contrange=FALSE
	}
	
	tfs = NULL 
	ntfs = NULL
	times = NULL
	ntimes = NULL
		
	if (contrange==TRUE){ # if all names have a TF_time structure
					
		if (length(grep(".+_.+", nms))==length(nms)){
						
			tfs = gsub("(.+)_.+","\\1",nms)
			tfs = unique(tfs)
			ntfs = length(tfs)

			times = gsub(".+_(.+)","\\1",nms) 
			times = unique(times)
			ntimes = length(times)

			if (ntimes*ntfs != length(nms)){
				contrange = FALSE
			}
		
		}
		else{
			contrange = FALSE
		}
		
	}
	
	if (contrange == FALSE){
				
		tfs = nms
		ntfs = length(tfs)
		times = ""
		ntimes = 1
	
	}
			
	if (!is.factor(cqbm[,classcol])){
		cqbm[,classcol] = as.factor(cqbm[,classcol])
	}
	else if (refactor==TRUE){
		
		if (is.integer(levels(cqbm[,classcol]))){
			cqbm[,classcol] = as.factor(as.integer(cqbm[,classcol]))
		}
		else{
			cqbm[,classcol] = as.factor(as.character(cqbm[,classcol]))
		}
	
	}
	#this refactoring is done because there may be unused levels 
	#in the input cqbm, as shown by experience
	
	rows = ceiling(length(levels(cqbm[,classcol]))/per.row)
	par(mfrow=c(rows,per.row))
	
					
	for(cc in levels(cqbm[,classcol])){
			
		m = matrix(data=numeric(0),nrow=ntimes,ncol=ntfs)
		ci.u = matrix(data=numeric(0),nrow=ntimes,ncol=ntfs)
		ci.l = matrix(data=numeric(0),nrow=ntimes,ncol=ntfs)
		mmin = matrix(data=numeric(0),nrow=ntimes,ncol=ntfs)
		mmax = matrix(data=numeric(0),nrow=ntimes,ncol=ntfs)
				
		for (j in 1:ntfs){
			for (i in 1:ntimes){
				if (ntimes>1){
					thissignal = cqbm[cqbm[,classcol]==cc,tfstcol+(j-1)*ntimes+i-1]
				}
				else{
					thissignal = cqbm[cqbm[,classcol]==cc,tfrange[j]]
				}
				        
				m[i,j] = mean(thissignal)
				if (length(thissignal)>1){	
					sdev = sd(thissignal)
				}
				else{
					sdev = 0
				}
				
				if (upperonly==T){
					if (m[i,j]>=0){
						ci.u[i,j] = m[i,j]+sdev
						ci.l[i,j] = m[i,j]
					}	
					else{
						ci.u[i,j] = m[i,j]
						ci.l[i,j] = m[i,j]-sdev
					}
				}
				else{
					ci.u[i,j] = m[i,j]+sdev
					ci.l[i,j] = m[i,j]-sdev	
				}
				
#				cat("\t",names(cqbm)[tfstcol+(j-1)*ntimes+i-1],"\n")
			}	
			
		}

		res[[as.character(cc)]] = m	
		res[[paste(as.character(cc),"_sd",sep="")]] = sdev
		res[[paste(as.character(cc),"_ci.u",sep="")]] = ci.u
		res[[paste(as.character(cc),"_ci.l",sep="")]] = ci.l
		
	}

	mmin = +Inf
	mmax = -Inf	

	if (commonylim == TRUE){
		for(cc in levels(cqbm[,classcol])){
			m = res[[as.character(cc)]]	
			sdev = res[[paste(as.character(cc),"_sd",sep="")]]
			
			thismin = ifelse(min(m)>0, 0, min(m+sdev))
			if (thismin<mmin){
				mmin = thismin
			}
			thismax = ifelse(max(m)>0, max(m+sdev), 0)
			if (thismax>mmax){
				mmax = thismax
			}	
		}
							
		mmin =  mmin*1.1
		mmax =  mmax*1.1

	}

	for(cc in levels(cqbm[,classcol])){
		main.bp = character(0)
		if (single==T){
			main.bp = cc
		}
		else{
			main.bp = paste("class",cc,"( n =",length(cqbm[cqbm[,classcol]==cc,1]),")")
		}

		m = res[[as.character(cc)]]	
		ci.u = res[[paste(as.character(cc),"_ci.u",sep="")]]
		ci.l = res[[paste(as.character(cc),"_ci.l",sep="")]]

		if (is.null(ylim)){
			ylim = c(mmin,mmax)
		}

		if (commonylim == TRUE){
			barplot2(m, beside=T, legend = times, plot.ci=plot.ci, ci.u = ci.u, ci.l = ci.l, names.arg = tfs, legend.text=legend.text,main = main.bp, ylim=ylim, ... ) 
		}
		else{
			barplot2(m, beside=T, legend = times, plot.ci=plot.ci, ci.u = ci.u, ci.l = ci.l, names.arg = tfs, legend.text=legend.text,main = main.bp, ... )
		}
	}

	invisible(res)
	
}


symmetrise = function(qq, tfrange, fill.with = 0){
	
	left=min(tfrange)-1
	right=length(qq[1,])-max(tfrange)

	init.colunmn.number = length(qq[1,])
	
	nms = names(qq)[tfrange]
	tfs = gsub("(.+)_.+","\\1",nms)
	tfs = unique(tfs)
	
	times = gsub(".+_(.+)","\\1",nms) 
	times = unique(times)

	firsttimes = sort(as.numeric(gsub("^(\\d+)\\.\\d+","\\1", times, perl=T)))
	
	for (tf in tfs){
		for (time in times){
			if (! paste(tf,time,sep="_") %in% names(qq)){
					qq=cbind(qq,rep(fill.with, length(qq[,1])))
					names(qq)[length(qq[1,])] = paste(tf,time,sep="_")
			}
		}
	}
	
#	print(head(qq))
	
	ordr = numeric(0)	
	if (left){ordr = c(1:left)}
	for(tf in tfs){
		for (firsttime in firsttimes){
			ordr = c(ordr, grep(paste("^",tf,"_",firsttime,"\\.",sep=""),names(qq),perl=T))
		}
	}
	if (right){ordr = c(ordr, (init.colunmn.number-right+1):init.colunmn.number)}
	
	qq <- qq[,ordr]
	
	print (head (qq))
	
	invisible(qq)
}


#deprecated synonims

cqbmReadPlotSave = function (...){ 
	chReadPlotSave(...)	
}

cbReadPlotOne = function (...){
	cbReadPlotClassSave(...)
}

cbPlotEach = function(...){
	cbPlotEachSave(...)
}

