## (c) Mikhail Spivakov 2014

library(gplots)
library(gtools)
library(RColorBrewer)
 
cluster_CRMs = function (cqbm, hmcol=colorRampPalette(c("#0000fe","black","#fdfd00"), space="Lab")(256), 
                         mycol = sample(rainbow(256)), key=F, density.info="none", 
                         tfrange=4:18, clustexcludenames=NULL, breaks=seq(-1,4,length.out=257), # tfrange=4:18 I have changed this
                         margins= c(12,3), dist = c("spearman", "euclidian")[2], 
                         clustmethod=c("complete","ward.D")[2], file="", tree.k = 10,  
                         clustcols=T,  rowpanel=T, rowdendr=T,coldendr=T, labRow="", ...){
 
 
 	thiscl = as.matrix(cqbm[ ,tfrange])
	
 	if (length(clustexcludenames)>0){
		whichexclude = which (attr(thiscl,"dimnames")[[2]]==clustexcludenames)
		thiscl2 = thiscl[, -whichexclude]
	}
	else{
		thiscl2 = thiscl
	}
	
	if(dist =="spearman"){
		dr = as.dist(1-cor(t(thiscl2),method="spearman"))
	}
	else{
		dr = dist(thiscl2)
	}
	
	hr = hclust(dr, method=clustmethod)

	if (rowpanel==T){
		mycl = cutree(hr, k=tree.k)  ## here
		mycol = mycol[as.vector(mycl)]
	}
	
	if (clustcols==T){
		if(dist =="spearman"){
			dr = as.dist(1-cor(thiscl,method="spearman"))
		}
		else{
			dr = dist(t(thiscl))
		}
		
		hc = hclust (dr,method=clustmethod)
		
		Colv = as.dendrogram(hc)
		if (rowdendr==T & coldendr==T){dendrogram = "both"}
		else {  if (rowdendr==T){dendrogram="row" }
				else  { 
					if (coldendr==T){ dendrogram="column"}
				}  
		}
	}
	else{
		Colv = FALSE
		if (rowdendr==T){dendrogram="row"}
		else{ dendrogram="none" }
	}


    if (rowpanel==T){
		cat("Row labels order (from bottom to top):", unique(mycl[labels(as.dendrogram(hr))]), "\n")
	    heatmap.3(thiscl, dendrogram=dendrogram, Colv=Colv, Rowv=as.dendrogram(hr), trace="none", labRow=labRow, RowSideColors=mycol, col=hmcol, key=key, density.info=density.info, breaks=breaks, margins=margins, rowclnames= rev(unique(mycl[labels(as.dendrogram(hr))])), ...)
    }
    else{
    	heatmap.3(thiscl, dendrogram=dendrogram, Colv=Colv, Rowv=as.dendrogram(hr), trace="none", labRow=labRow, col=hmcol, key=key, density.info=density.info, breaks=breaks, margins=margins, rowclnames= rev(unique(mycl[labels(as.dendrogram(hr))])), ...)
    }    
    
    if (file!=""){
    	readline("Choose a nice window shape and press any key...")
		quartz.save(file, type="pdf")
	}
    
   	if (rowpanel==T){
   		invisible(cbind(cqbm, mycl))
   	}
}


plot_this_CRM_heatmap = function(cqbm, classOrder, CexAngleColourPerClass, clcol=which(names(cqbm)=="V3"), 
							legendMar=13, ...){
	
	cfqbm = cqbm
#	cfqbm$V3 = factor(cqbm[,clcol], levels=classOrder)

	CexAngleColourPerClass[,1] = factor(CexAngleColourPerClass[,1], levels = classOrder)
	CexAngleColourPerClass = CexAngleColourPerClass[order(CexAngleColourPerClass[,1]),  ]
	plot_CRM_heatmap(cfqbm, classcol=clcol, classOrder=classOrder,
							rowSideTextCex = CexAngleColourPerClass[,2],  
							rowSideTextAngle = CexAngleColourPerClass[,3], 
							mycol = CexAngleColourPerClass[,4], 
							legendMar=legendMar, ...)			
}

#THIS FUNCTION
plot_CRM_heatmap = function(cqbm, hmcol=colorRampPalette(c("#0000fe","black","#fdfd00"), space="Lab")(256), classcol = which(names(cqbm)=="V3"),
                            mycol = sample(rainbow(256)), key=TRUE, density.info="none", tfrange=7:16, clustexcludenames=NULL,
							excludeclasses=NULL, breaks=seq(-1,4,length.out=257), 
							tmax=NULL, tmin=NULL, 
							# can also impose limits just using breaks, unless the clustering doesn't look nice this way
       						margins= c(12,3), dist = c("spearman", "euclidian")[2], clustmethod=c("complete","ward.D")[2], file="", 
       						clustercols=FALSE, coldist = c("spearman", "euclidian")[2], colclustmethod=c("complete","ward.D")[2],
       						separateRows="RowSide", seplcol= "white", writeClassNames = TRUE, classOrder=NULL, ...){
		
	mycl = NULL
	allcl = matrix(ncol=length(tfrange), nrow=0)
	
	
	classes = NULL
	
#	if(all(is.na(as.numeric(as.character(cqbm[,classcol])))==FALSE)){
#		cat("all class levels are numeric\n")
#		classes = sort(unique(as.numeric(as.character(cqbm[,classcol]))))
#	}
#	else{
#		cat("non-numeric class names found\n")
#		classes = levels( as.factor(cqbm[, classcol]) )
#	}
	
	## excludeclasses <- c(11,15)
	cqbm = cqbm[ !cqbm[,classcol]%in%excludeclasses , ] ##removing any unwanted class/cluster from the table of interest 
	
	if (missing(classOrder)){ #by default yes
		classOrder = unique(cqbm[,classcol]) #order of names of clusters appearing for first time in table of interest
	}
	
	cqbm[, classcol] = factor(cqbm[, classcol], levels=classOrder) #as factor the column telling us the cluster the interaction pertains to
	
	cat("total len: ",length(cqbm[,classcol]),"\n")

	#print(unique(cqbm[,classcol]))
	
	classes = levels(cqbm[, classcol]) #what clusters exist? #what names do they have? 

	ColV = FALSE
	dgr = "none"
	
	if (clustercols==TRUE){ #not by default
		if (coldist=="spearman"){
			dc = as.dist(1-cor(cqbm[,tfrange], method="spearman"))
		}
		else{
			dc = dist(t(cqbm[tfrange]))
		}
		hc = hclust(dc, method=colclustmethod)
		ColV = as.dendrogram(hc)
		dgr = "column"
	}
	

	if (!is.null(tmax)){ #not by default
		cqbm[,tfrange] = apply(cqbm[,tfrange],c(1,2),function(x)if(!is.na(x)){if(x<tmax){return(x)}else{return(tmax+1)}}else{return(tmax+1)})
	}
	if (!is.null(tmin)){ #not by default
		cqbm[,tfrange] = apply(cqbm[,tfrange],c(1,2),function(x)if(!is.na(x)){if(x<tmin){return(x)}else{return(tmin-1)}}else{return(tmin-1)})
	}


	i=1
	#classes <- classes[c(1:16)]  ME Andrea --- removing those clusters with only one interaction
	for (cl in classes ){ #classes = levels(cqbm[, classcol]) ; for each cluster: 
				
		cat("Processing class = ", as.character(cl),"\n")
		
		thiscl = as.matrix(
					cqbm[as.factor(cqbm[,classcol])==cl, #extracts cluster of interactions of interest
					tfrange]
				 )
								
		if (length(clustexcludenames)>0){ 
			whichexclude = which (attr(thiscl,"dimnames")[[2]]==clustexcludenames)
			thiscl2 = thiscl[, -whichexclude]
		}
		else{
			thiscl2 = thiscl
		}
		
		if(dist =="spearman"){
			dr = as.dist(1-cor(t(thiscl2),method="spearman"))

		}
		else{
			dr = dist(thiscl2)
		}
		hr = hclust(dr, method=clustmethod) #clusters interactions
		rowInd = order.dendrogram(as.dendrogram(hr)) #obtains order of interactions
				
		thiscl = thiscl[rowInd,] # reorders interactions (rows)
		
		cat("\tlength thiscl=", length(thiscl[,1]),"\n")		 	
		
		allcl = rbind(allcl, thiscl) #allcl in the end is the same matrix as input qbmc but with reordered rows by clustering of clusters 

		mycl = c(mycl, rep(i, length(thiscl[,1]))) # create a 1..n index system for classes instead of bothering with factor levels
		i=i+1
	}
	
	cat("Length allcl=",length(allcl[,1]),"\n")
	
    mycol = mycol[as.vector(mycl)]        
    heatmap.3(allcl, dendrogram=dgr, Colv=ColV, Rowv=F, trace="none", labRow="", RowSideColors=mycol, col=hmcol, key=key, density.info=density.info, breaks=breaks, margins=margins, rowclnames=classes, separateRows=separateRows, seplcol=seplcol, writeClassNames = writeClassNames, ...)

    if (file!=""){
      #quartz or sth!!
	}
}


heatmap.3 = function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, col = "heat.colors", colsep, 
    rowsep, sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, 
    notecex = 1, notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5,  denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE), densadj = 0.25, main = NULL, xlab = NULL, 
    ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
    density.info = c("histogram", "density", "row.legend", "none"), # row.legend added in heatmap.3
    # the parameters below are new in heatmap.3
    legendMar = NULL, legendTextCex = 1, legendText = NULL,
    rowclnames=NULL, writeClassNames=FALSE, rowSideTextAngle = 0, rowSideTextCex = 1,
    separateRows=NULL, seplwd=2, seplcol="red", ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if ((Colv == "Rowv") && (!isTRUE(Rowv) || is.null(Rowv))) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) 
        if (missing(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    if (length(breaks) == 1) {
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
            length = breaks)
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    else if (is.character(col) && length(col) == 1) 
        col <- do.call(col, list(ncol))
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[] <- ifelse(x < min.breaks, min.breaks, x)
    x[] <- ifelse(x > max.breaks, max.breaks, x)
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    #layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)

		# new in heatmap.3
		if (writeClassNames==TRUE){
			 if (length(rowclnames)==length(unique(RowSideColors))){
				prevClassesLen = 0 
		 	    i = 1
			 	for (colour in unique(RowSideColors)){
					
					rss = numeric(0) 
					
					if (length(rowSideTextAngle)==1) rss = rowSideTextAngle 
					else rss = rowSideTextAngle[i] 
	
					rsc = numeric(0)

					if (length(rowSideTextCex)==1) rsc = rowSideTextCex 
					else rsc = rowSideTextCex[i] 
				
			 		thisClassLen = length(RowSideColors[RowSideColors==colour])
			 				 		
        			text(0,(length(rowInd)-prevClassesLen-thisClassLen/2)/length(rowInd), rowclnames[i], srt=rss, cex=rsc)
        			prevClassesLen = prevClassesLen+thisClassLen
			 		i = i+1
			 	}
			 }
			 else{
			 	stop ("rowclnames should have the same length as the number of classes.\ncheck if mycol setting is correct")
			 } 
		}
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none") {
        x <- t(x)
        cellnote <- t(cellnote)
    }
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    
    
  #  print (length(breaks))
#    print (length(col))
        
    ##################### main image ###########################
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    ############################################################
    
    if (!gtools::invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
            length(csep)), xright = csep + 0.5 + sepwidth[1], 
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        for (i in colInd) {
            if (!is.null(vline)) {
                vline.vals <- scale01(vline, min.scale, max.scale)
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        for (i in rowInd) {
            if (!is.null(hline)) {
                hline.vals <- scale01(hline, min.scale, max.scale)
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    
    # Added in heatmap.3
     
    if (!missing(separateRows)){
    	if (length (separateRows)==1 && "RowSide" %in% separateRows ){
    		i = 0
    		for (colour in unique (RowSideColors)){
    			i = i + length (RowSideColors[RowSideColors==colour])
    			#cat ("i=",i,"\n")
	    		abline(h=length(rowInd)-i, lty=1, lwd=seplwd, col=seplcol) # separateRows numbers rows from the top!
    		}
    	}
    	else for (i in separateRows){
    		abline(h=length(rowInd)-i, lty=1, lwd=seplwd, col=seplcol) # separateRows numbers rows from the top!
    	}
    }
    
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])

    if (key) {
        if (density.info!="row.legend"){
			    
			    # new in heatmap.3
	   			if (missing(legendMar)){
				    if (! density.info %in% c("density","histogram"))  legendMar = 4
				  	else legendMar = 2
			  	}
		    	par(mar = c(5, 4, legendMar, 1), cex = 0.75) # in heatmap.2 : 5,4,2,1 in all cases
			  	
			  	
      	       if (symkey) {
	            	max.raw <- max(abs(x), na.rm = TRUE)
         		   min.raw <- -max.raw
        		}
        		else {
            		min.raw <- min(x, na.rm = TRUE)
            		max.raw <- max(x, na.rm = TRUE)
        		}
        		z <- seq(min.raw, max.raw, length = length(col))
        		image(z = matrix(z, ncol = 1), col = col, breaks = breaks, 
            		xaxt = "n", yaxt = "n")
        		par(usr = c(0, 1, 0, 1))
        		lv <- pretty(breaks)
        		xv <- scale01(as.numeric(lv), min.raw, max.raw)
        		axis(1, at = xv, labels = lv)
        		if (missing(legendText)){     	    
					if (scale == "row") 
        				legendText = "Row Z-Score"
	        		else if (scale == "column") 
		    	    	legendText = "Column Z-Score"
        			else legendText = "Value"
        		}
        		
           		mtext(side = 1, legendText, line = 2, cex=legendTextCex)
        }

        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        
        # Added in heatmap.3
        else if (density.info == "row.legend"){
     	 legend("center", legend = rowclnames,
		        col = rev(unique (RowSideColors[rowInd] )), pch=15, cex=0.8, horiz=T)
#		        pt.cex=2,
#		        lty = c(0, 0, 2, 2))
			  cat("ColSide colours from top to bottom:\n")
        	  print (rev(unique(RowSideColors[rowInd])))
        	}
 
#        else if (notitle)title("\nColor Key") removed in heatmap.3
    }
    
    else plot.new()
    invisible(list(rowInd = rowInd, colInd = colInd))
}

