## (c) Mikhail Spivakov 2011

source("/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/plot_class_heatmap.R")

templates = "/gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/programs/autoclass/autoclass/templates/"

autoclass = function (data,
			scale=TRUE,
			delog=FALSE,
			zero_point = "min", 
			rel_error = 0.1,  
			# for type "real scalar"; a single value or one for _each_ _column_! 
			# (only those corresponding to columns of type real scalar will be used)  
			error = 1, 
			# for type "real location"; a single value or one for _each_ _column_! 
			# (only those corresponding to columns of type real scalar will be used)
			single_normal_cn=NULL, # column numbers in data
			single_multinomial=NULL, # column numbers in data
			single_normal_cm=NULL, # column numbers in data
			multi_normal_cn=NULL, # column numbers in data
			id=NULL, # for autoclass2tab.pl id:file=<db2-fname>,col=<id>. If not provided, will be created automatically
			# will automatically set all the non-modelled columns to the type ignore
			dir = paste(getwd(),"/", project, sep=""),
			project="sample",
			types = "real scalar", # will automatically substitute this with "dummy nil" for id and ignore
			runautoclass=TRUE, 
			runautoclass2tab=TRUE,
			pmin = 0.5, # NULL if should not apply
			rmin = 2, # NULL if should not apply
			plot=TRUE, 
			plotrange="includedcolumns", # which columns to plot. Default = all includedcolumns
						   # if names have the structure TF_time, should have an equal amount of times per column
						   # use symmetrise() from plot_class_heatmap.R if this is not the case
			interactive=TRUE, # should the user be asked to choose a nice window shape
			starts =  c(2, 3, 5, 7, 10, 15, 25),
			fn ="converge_search_3",
			max_cycles=1000,
			max_n_tries=100, ...){ 	# ... - arguments to pass to plot_class_barcharts and through there to barplot2

	# some last-minute aliases
	output_dir = dir
	project_name = project
	start_j_list = starts
	try_fn_type = fn
	
	datavarname = as.character(substitute(data))
	
	ignore = NULL
	
	optwarn = options("warn")$warn
	options(warn=-1)
		
	if (missing (single_normal_cn) & missing(single_multinomial) & missing(single_normal_cm) & missing(multi_normal_cn)){
		stop("No model specified for the data. Attributes can be single_normal_cn, single_multinomial, single_normal_cm or multi_normal_cn")
	}

	cat("\nPreparing data...\n\n")

	if (!is.data.frame(data)){
		data = as.data.frame(data)
	}
	
	ncol = length(data[1,])

	initid = id # only needed to correctly advise what command to execute if want to use the 2nd-best classification
	# if Id is missing, add an index column
	if (missing(id)){
		cat("Adding an id column...\n")
		data[, "autoclassid"] = 1:length(data[,1])
		ncol = ncol+1
		id = ncol	
		initid = "NULL"
	}
	else{
		cat("The id column is set to", id, "\n")
	}
	
	# Save the dataframe before exponentiating/scaling, but after adding the ID column if missing
	original.data = data
		
	includedcolumns = sort(c(single_normal_cn, single_multinomial, single_normal_cm, multi_normal_cn)) 
	if (plotrange == "includedcolumns"){ # ie, if the range of columns is continuous
		plotrange = includedcolumns
	}
	
	whichidentical = sapply(includedcolumns, function(x)sd(data[,x]))

	remove = NULL
	for (i in 1:length(whichidentical)){
		if (! whichidentical[i]){
			cat ("Column", includedcolumns[i], "has a stdev of zero and will be ignored\n")
			remove = c(remove, includedcolumns[i])
		}
	} 
	if (length(remove)>0){
		includedcolumns = includedcolumns[-which(includedcolumns %in% remove)]
		single_normal_cn = single_normal_cn[-which(single_normal_cn %in% remove)]
		single_multinomial = single_multinomial[-which(single_multinomial %in% remove)]
		single_normal_cm = single_normal_cm[-which(single_normal_cm %in% remove)]
		multi_normal_cn = multi_normal_cn[-which(multi_normal_cn %in% remove)]
	}

	if ( ! all(includedcolumns == 1:length(data[1,]))){
		ignore = (1:length(data[1,]))[! (1:length(data[1,])) %in%  includedcolumns  ]
		cat ("The following columns will be ignored in clustering:", paste(ignore, collapse=", "), "\n")
	}


	if (any(apply(data[, single_normal_cm], 2, function(x)any(is.na(x)))==TRUE)){
		stop("Missing values found for columns whose model type is not signle_normal_cm")
	}
		

	if (! all(types%in%c("real scalar","real location", "discrete nominal", "dummy nil"))){
		stop ("Unsupported attribute types specified. Supported types are: \"real scalar\", \"real location\", \"discrete nominal\", \"dummy nil\"")
	}
	
	if (! try_fn_type %in% c("converge_search_3","converge_search_4","converge")){
		stop ("Unsupported fn type (try_fn_type). Supported types: \"converge_search_3\",\"converge_search_4\",\"converge\"")
	}
	

	if (length(types)==1 & ncol!=1){
		types = rep(types,ncol)
	}	
	if (length(ignore)>0){
		types[ignore] = "dummy nil"
	}

	real.scalar = which(types=="real scalar")
	if (delog == TRUE){
		for (i in real.scalar){
			data[,i] = exp(data[,i])
		}
	}

  if (scale == TRUE){
		for (i in real.scalar){
			data[,i] = scale(data[,i])
		}
	}
  
	for (i in single_normal_cm){
	  data[is.na(data[,i]),i] = "?"
	}
	
	# The below will generate a warning for non-numeric NAs 
	if (length(grep("min", zero_point))>0){ 
		mins = sapply(1:ncol,function(i)if (is.numeric(data[,i]) | is.integer(data[,i])){min(data[,i])}else{NA})
		mins[is.na(mins)] = 0 # a patch for ignored columns of type character
	}
	options(warn=optwarn)

	if (length(zero_point)==1 & ncol!=1){
		zero_point = rep(zero_point, ncol)
	}
	for (i in 1:ncol){
		if(zero_point[i]=="min"){
			zero_point[i] = mins[i]
		}
	}
	if (length(zero_point)<length(types[types=="real scalar"])){
		stop ("Zero point should be specified for each or all attributes of type real scalar")
	}
		
	if (length(rel_error)==1 & ncol!=1){
		rel_error = rep(rel_error, ncol)
	}
	if (length(rel_error)<length(types[types=="real scalar"])){
		stop ("Relative error should be specified for each or all attributes of type real scalar")
	}
	
	if (length(error)==1 & ncol!=1){
		error = rep(error, ncol)
	}
	if (length(error)<length(types[types=="real location"])){
		stop ("Relative error should be specified for each or all attributes of type real location")
	}
	
	isFactor = sapply(1:ncol,function(x)is.factor(data[,x]))
	whichfactors = which(isFactor==TRUE)
	whichmeaningfulfactors = whichfactors[! whichfactors %in% ignore]
	if(length(whichmeaningfulfactors)>0){
		if (any(types[whichmeaningfulfactors] != "discrete nominal") | any(whichmeaningfulfactors != single_multinomial)){
			stop("For discrete factors that should not be ignored, only the discrete nominal data type is supported and only the single multinomial data model")
		}	
		for (i in whichmeaningfulfactors){
			data[,i] = as.integer(data[,i])
		}
		cat ("Factor column(s)", paste(whichmeaningfulfactors, collapse=", "), "were converted to integers according to factor levels\n")
	}
	
	cat ("\nWriting autoclass input files...\n")
	
	# Create project folder. Warn if exists

	if (file.exists(output_dir)){
		go = readline(paste("Folder", output_dir, "exists. Confirm overwrite [y/N]"))
		if (go == "y"){
			unlink(output_dir, recursive=TRUE)
		}
		else{
			stop ("Folder overwrite cancelled")
		}
	}
	ok = dir.create(output_dir)
	if (ok==FALSE){
		stop("mkdir returned an error")
	}

	# Change into the project folder
	 
	curr.dir <- getwd()
	setwd(output_dir)

	# Copy the template of .db2 file (renamed)
	
	db2.name = paste(project_name,".db2", sep="")
	ok = file.copy(paste(templates,"sample.db2",sep=""), db2.name)
	if(ok==FALSE){
		stop ("copying of sample.db2 to", paste(getwd(), db2.name,sep="/"), "returned an error")
	}
	
	# Append data to the .db2 file
	write.table(data, file=db2.name, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

	# Copy the template of .s-params file (renamed)
	s.params.name = paste(project_name,".s-params", sep="")
	ok = file.copy(paste(templates,"sample.s-params",sep=""), s.params.name)
	if(ok==FALSE){
		stop ("copying of sample.s-params to", paste(getwd(), s.params.name,sep="/"), "returned an error")
	}

	# Append to the .s-params file:
	
	s.params = file(s.params.name, "a")
	cat("\ntry_fn_type = \"", try_fn_type, "\"\n", file=s.params, sep="")
	cat("max_cycles =", max_cycles, "\n", file=s.params)
	cat("max_n_tries =", max_n_tries, "\n", file=s.params)
	cat("start_j_list =", paste(start_j_list, collapse = ", "), "\n", file=s.params)
	close(s.params)
	
	# Copy the template of the .hd2 file (renamed)

	hd2.name = paste(project_name,".hd2", sep="")
	ok = file.copy(paste(templates,"sample.hd2",sep=""), hd2.name)
	if(ok==FALSE){
		stop ("copying of sample.hd2 to", paste(getwd(), hd2.name,sep="/"), "returned an error")
	}

	# Append the attribute descriptions to the .hd2 file
	hd2 = file(hd2.name, "a")
	cat("\n\n", file=hd2)
	cat("number_of_attributes", ncol, "\n\n", file=hd2)
	for (i in 1:ncol){
		if (types[i]=="dummy nil"){
			cat(i-1, " dummy nil \"", names(data)[i], "\"\n", sep="", file=hd2)
		}
		else if (types[i]=="real scalar"){
			cat(i-1, " real scalar \"", names(data)[i], "\" zero_point ", zero_point[i], " rel_error ", rel_error[i], "\n", sep="", file=hd2)
		}
		else if (types[i]=="real location"){
			cat(i-1, " real location \"", names(data)[i], "\" error ", error[i], "\n", sep="", file=hd2)
		}
		else if (type[i]=="discrete nominal"){
			cat(i-1, " discrete nominal \"", names(data)[i], "\" range ", max(data[,i]), "\n", sep="", file=hd2)
		}
	}
	close(hd2)
	
	# Copy the template of the .model file (renamed)
	
	model.name = paste(project_name,".model", sep="")
	ok = file.copy(paste(templates,"sample.model",sep=""), model.name)
	if(ok==FALSE){
		stop ("copying of sample.model to", paste(getwd(), model.name,sep="/"), "returned an error")
	}
	
	# Append the model descriptions to the .model file 
	model = file(model.name, "a")	
	nmodels = 0
	model.string=""
	if (length(ignore)>0){
		model.string=paste(model.string, "ignore ",  paste(ignore-1, collapse=" "), "\n", sep="")
		nmodels = nmodels+1
	}
  if (length(single_normal_cn)>0){
		model.string=paste(model.string, "single_normal_cn ", paste(single_normal_cn-1, collapse=" "), "\n", sep="")
		nmodels = nmodels+1
	}
	if (length(single_normal_cm)>0){
		model.string=paste(model.string, "single_normal_cn ", paste(single_normal_cn-1, collapse=" "), "\n", sep="")
		nmodels = nmodels+1	
	}
	if (length(single_multinomial)>0){
		model.string=paste(model.string, "single_multinomial ", paste(single_multinomial-1, collapse=" "), "\n", sep="")
		nmodels = nmodels+1	
	}
	if (length(multi_normal_cn)>0){
		model.string=paste(model.string, "multi_normal_cn ", paste(multi_normal_cn-1, collapse=" "), "\n", sep="")
		nmodels = nmodels+1
	}
		
	cat("\n", file=model)
	cat ("model_index", 0, nmodels, "\n\n", file=model)
	cat (model.string, file=model)
	close (model)
	
		
	# Copy the .r-params file (renamed)
	
	r.params.name = paste(project_name,".r-params", sep="")
	ok = file.copy(paste(templates,"sample.r-params",sep=""), r.params.name)
	if(ok==FALSE){
		stop ("copying of sample.r-params to", paste(getwd(), r.params.name,sep="/"), "returned an error")
	}
	
	if (runautoclass==FALSE){
		cat("Runautoclass set to FALSE, exitting. You can check the files prepared for autoclass in", output_dir, "\n")
		return(0)
	}
	
	# Run autoclass in the -search mode
	cat ("\nRunning autoclass in the -search mode. This might take a while and little output may appear straight away...\n\n")
	command = paste(autoclass.path, "-search", db2.name, hd2.name, model.name, s.params.name)
	cat (command, "\n")
	err = system(command)
	if (err){
		stop ("Autoclass -search returned an error")
	}	

	# Run autoclass in the -report mode
	cat ("\nRunning autoclass in the -report mode...\n\n")
	results.name = paste(project_name,".results-bin", sep="")
	search.name = paste(project_name,".search", sep="")
	command = paste(autoclass.path, "-reports", results.name, search.name, r.params.name)
	cat (command, "\n")
	err = system(command)
	if (err){
		stop ("Autoclass -report returned an error")
	}	

	# Unless runautoclass2tab=FALSE, run autoclass2tab for the two variants and plot the results
	if (runautoclass2tab==TRUE){
		cat("\nAnalysing the best match...\n")
		cat("autoclass2tab(data=", datavarname, ", project=\"", project_name,"\", pmin=", pmin, ", rmin=", rmin, ", id=", id,", whichres = 1, plot=", plot, ", plotrange = ", deparse(plotrange), ", interactive=",interactive,", ...)\n", sep="")
		
		cqbm = autoclass2tab(data=original.data, project=project_name, pmin=pmin, rmin=rmin, id=id, whichres = 1, plot=plot, plotrange=plotrange, interactive=interactive, ...) # save cqbm for the first best classification
		
		cat("\nAnalysing second-best match...\n")
		cat("autoclass2tab(data=", datavarname, ", project=\"", project_name,"\", pmin=", pmin, ", rmin=", rmin, ", id=", id,", whichres = 2, plot=", plot, ", plotrange = ", deparse(plotrange), ", interactive=",interactive,", ...)\n", sep="")
		
		autoclass2tab(data=original.data, project=project_name, pmin=pmin, rmin=rmin, id=id, whichres = 2, plot=plot, plotrange=plotrange, interactive=interactive, dir=dir, ...) # only draw the second best classification
	}

	if(! is.null(curr.dir)){
		setwd(curr.dir)
	}
	
	if (runautoclass2tab==TRUE){
		cat("\nReturning the best match.\nFor the second-best match, save the result of:\nautoclass2tab( ", datavarname, ", project=\"", project_name, "\", pmin=", pmin, ", rmin=", rmin, ", id=", initid, ", whichres=2", ", dir=\"", output_dir,  "\", plot=FALSE)\n\n", sep="")
		invisible (cqbm)
	}
}

autoclass2tab = function(data=data, project="sample", pmin=0.5, rmin=2, id=NULL, whichres = 1, dir = getwd(), plot = TRUE, plotrange=NULL, interactive=TRUE, ...){ 
	# ... - arguments to pass to plot_class_barcharts and through there to barplot2
 	oldd = getwd()
 	project_name = project
	setwd(dir)
	
  data =as.data.frame(data) # in case it was a matrix
    
	print (head(data))
	#print(plotrange)
	
	cat("Running autoclass2tab for classification", whichres, "....\n")
	
	report.class.data.name = paste(project_name, ".class-data-", whichres, sep="")
	command = paste (autoclass2tab.path, report.class.data.name)
	if (!is.null(pmin)){
		command = paste(command, " pmin=", pmin, sep="")
	}
	if (!is.null(rmin)){
		command = paste(command, " rmin=", rmin, sep="")
	}
	if (!is.null(id)){
		command = paste(command, " id:file=", project_name, ".db2,col=", id, sep="")
	}
	else{
		cat ("adding ID as last column- this assumes that it was added as well when running autoclass...")
		data$autoclassid = 1:length(data[,1])
		id = which(names(data)=="autoclassid")	
	}
	tab.name = paste(project_name, "-", whichres, ".tab", sep="")

	print (tab.name)
	command = paste (command, ">", tab.name, sep="")
 	
	cat (command, "\n")
	err = system (command)
	if (err){
		stop ("Autoclass2tab returned an error")
	}

  #classdata = make_cqbm(class_tab_fname = tab.name, qbm = data, qbm_idname = names(data)[id])
  #return(classdata)

	if (!is.null(plotrange) & plot==TRUE ){
		plot_class_barcharts(classdata, classcol=3, tfrange=which(names(classdata)%in%names(data)[plotrange]), ...) 

		plot.name = paste(project_name, "-", whichres, ".pdf", sep="") 
		if (interactive==TRUE){
			readline("Choose a nice window shape and press any key...\n")
		}
		# quartz.save(plot.name, "pdf")
	}

	setwd(oldd)
	#invisible(classdata)		
	
}
# 
