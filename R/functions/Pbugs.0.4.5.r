###############################################################
## Esta función se corresponde con la versión 0.4.5 de Pbugs ##
###############################################################

require(lattice)
require("snow")
require("R2WinBUGS")
#require("R2jags")
 
#setClass("bugs")
#setClass("WinBUGS", contains="bugs")
#setClass("OpenBUGS", contains="bugs")
plot.OpenBUGS<-function(x, ...){R2OpenBUGS:::plot.bugs(x, ...)}
print.OpenBUGS<-function(x, ...){R2OpenBUGS:::print.bugs(x, ...)}
plot.WinBUGS<-function(x, ...){R2WinBUGS:::plot.bugs(x, ...)}
print.WinBUGS<-function(x, ...){R2WinBUGS:::print.bugs(x, ...)}
 
Pbugs<-function(program = c("winbugs", "openbugs"), cluster=NULL, pbugs.directory=ifelse(.Platform$OS.type == "windows", "c:/Pbugs", path.expand("~/.wine/dosdevices/c:/Pbugs")), ...){
	program<-tolower(program)
	program <- match.arg(program)
	if(program=="winbugs"){
		require("R2WinBUGS")
		bugs.obj<-PWinBUGS(program="winbugs", cluster=cluster, pbugs.directory=pbugs.directory, ...)
#		bugs.obj<-PWinBUGS(program="winbugs", ...)
	}
	else{
		require("R2OpenBUGS")
		bugs.obj<-POpenBUGS(cluster=cluster, pbugs.directory=pbugs.directory, ...)
	}
	bugs.obj
}
 
 
plotHistory<-function (x,criteria=list(Rhat=1.1,n.eff=100),variables=NULL,mfrow=c(3,2),MaxPlots=10, ...){
	#if(!(class(x)%in%c("bugs","WinBUGS","OpenBUGS"))){stop(substitute(x)," must be an object of one of the classes: bugs, WinBUGS or OpenBUGS")}
  if(!("bugs"%in%class(x))){stop(substitute(x)," must be an object of one of class bugs")}
	if(is.null(variables)){
		SelVar<-which(x$summary[,"Rhat"]>criteria$Rhat | x$summary[,"n.eff"]<criteria$n.eff)
		if(length(SelVar)==0){cat("All variables accomplish the required convergence criteria. No variables have been plotted.\n");return(invisible())}
	}
	else{
		if(!is.vector(variables) | !is.character(variables)) stop("Variables should be a vector with the expressions of the variables to be plotted")
		aux<-unique(sapply(variables,function(y){strsplit(y,"\\[")[[1]][1]})) %in% x$root.short
		if(!all(aux)) stop(paste("Expression(s) ",unique(sapply(variables,function(x){strsplit(x,"\\[")[[1]][1]}))[which(aux==FALSE)]," does not correspond to any variable(s) in ",substitute(x),"\n",sep=""))
		allnames<-character()
		for(i in 1:length(variables)){
			if(length(grep("\\[",variables[i]))==0){
				this<-eval(parse(text=paste("x$mean$",variables[i],sep="")))
				if(is.vector(this) & length(this)==1) allnames<-c(allnames,variables[i])
				if(is.array(this)){
					auxdim<-dim(this)
					if(length(auxdim)==1){
						allnames<-c(allnames,paste(variables[i],"[",1:auxdim,"]",sep=""))
					}
					if(length(auxdim)==2){
						allnames<-c(allnames,paste(variables[i],"[",rep(1:(auxdim[1]),each=auxdim[2]),",",rep(1:(auxdim[2]),auxdim[1]),"]",sep=""))
					}
					if(length(auxdim)==3){
						allnames<-c(allnames,paste(variables[i],"[",rep(1:(auxdim[1]),each=prod(auxdim[-1])),",",rep(rep(1:(auxdim[2]),auxdim[1]),each=auxdim[3]),",",rep(1:(auxdim[3]),prod(auxdim[-3])),"]",sep=""))
					}
					if(length(auxdim)==4){
						allnames<-c(allnames,paste(variables[i],"[",rep(1:(auxdim[1]),each=prod(auxdim[-1])),",",rep(rep(1:(auxdim[2]),auxdim[1]),each=prod(auxdim[3:4])),",",rep(rep(1:(auxdim[3]),prod(auxdim[1:2])),each=auxdim[4]),",",rep(1:(auxdim[4]),prod(auxdim[-4])),"]",sep=""))
					}
					if(length(auxdim)>4) stop("plot.History does not support arrays with 5 or more dimensions")
				}
			}
			else{
				auxdim<-dim(eval(parse(text=paste("x$mean$",strsplit(variables[i],"\\[")[[1]][1],sep=""))))
 
				if(length(auxdim)==1){
					subset<-strsplit(strsplit(variables[i],"\\[")[[1]][2],"\\]")[[1]]
					if(substr(subset,1,1)=="-"){subset<-(1:(auxdim[1]))[eval(parse(text=subset))]}
					else{subset<-eval(parse(text=subset))}
					allnames<-c(allnames,paste(strsplit(variables[i],"\\[")[[1]][1],"[",subset,"]",sep=""))
				}
				else{
					presubset<-strsplit(strsplit(variables[i],"\\[")[[1]][2],"\\]")[[1]]
					aux<-strsplit(presubset,",")[[1]]
					presubset2<-list()
 
					for(j in 1:length(auxdim)){
						presubset2[[j]]<-aux[1]
						aux<-aux[-1]
						if(presubset2[[j]]=="") presubset2[[j]]<-1:(auxdim[j])
						else{
							if(substr(presubset2[[j]],1,1)=="-") presubset2[[j]]<-(1:(auxdim[j]))[eval(parse(text=presubset2[[j]]))]
							else{
								while(length(grep("\\(",presubset2[[j]]))==1 & length(grep("\\)",presubset2[[j]]))==0){presubset2[[j]]<-paste(presubset2[[j]],aux[1],sep=",");aux<-aux[-1]}
								if(class(try(eval(parse(text=presubset2[[j]])),silent=T))=="try-error"){stop(paste("Expression not supported: ", variables[i]))}
								presubset2[[j]]<-eval(parse(text=presubset2[[j]]))
							}
						}
					}
 
					subset<-as.vector(t(outer(presubset2[[1]],presubset2[[2]],function(x,y){paste(x,y,sep=",")})))
					if(length(auxdim)==3 | length(auxdim)==4) subset<-as.vector(t(outer(subset,presubset2[[3]],function(x,y){paste(x,y,sep=",")})))
					if(length(auxdim)==4) subset<-as.vector(t(outer(subset,presubset2[[4]],function(x,y){paste(x,y,sep=",")})))
					if(length(auxdim)>4) stop("plot.History does not support arrays with 5 or more dimensions")
 
 					allnames<-c(allnames,paste(strsplit(variables[i],"\\[")[[1]][1],"[",subset,"]",sep=""))
				}
			}
		}
		SelVar<-match(allnames,dimnames(x$summary)[[1]])
		if(any(is.na(SelVar))) stop("Expression(s) ", allnames[which(is.na(SelVar))]," does not correspond to any variable recorded in ",substitute(x),"\n")
	}
 
	if(length(SelVar)>prod(mfrow)){cat("The number of graphs to be plotted excedes the number of panels in the active window.\nPress return to scroll forth to the next set of plots.\n")}
	Sel<-x$sims.array[,,SelVar]
	par.old<-par()$mfrow
	par(mfrow=mfrow)
	on.exit(par(mfrow=par.old))
	for(i in 1:(dim(Sel)[3])){
		if(i==(prod(mfrow)+1)){
			par.ask<- par()$ask
			par(ask=TRUE)
  	  on.exit(par(ask=par.ask))
		}
		plot(Sel[,1,i],type="l",ylim=c(min(Sel[,,i]),max(Sel[,,i])),xlab="Iteration",ylab="",...)
		for(j in 2:(dim(Sel)[2])){
			lines(Sel[,j,i],type="l",col=j)
		}
		title(dimnames(Sel)[[3]][i])
		if(i==MaxPlots*prod(mfrow) & (dim(Sel)[3]>i)){
			cat("Maximum amount of plots allowed has been reached, increase MaxPlots if you still want to watch further plots\n")
			return(invisible())
		}
	}		
}
 
 
 
PWinBUGS<-function(data, inits, parameters.to.save, model.file, n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter/2), n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/n.sims)), n.sims = 1000, bin = (n.iter - n.burnin)/n.thin, debug = FALSE, DIC = TRUE, digits = 5, codaPkg = FALSE, bugs.directory = ifelse(.Platform$OS.type == "windows", "c:/Program Files/WinBUGS14", path.expand("~/.wine/dosdevices/c:/Program Files/WinBUGS14")), program = c("WinBUGS", "OpenBUGS", "winbugs", "openbugs"), cluster=cluster, pbugs.directory=pbugs.directory, working.directory = NULL, clearWD = FALSE, useWINE = .Platform$OS.type != "windows",WINE="/usr/bin/wine",newWINE=TRUE, WINEPATH="/usr/bin/winepath", bugs.seed = NULL, summary.only = FALSE, save.history = !summary.only, over.relax = FALSE){
 
	#####################
	# Pbugs-specific code
	if(summary.only) {
		summary.only<-FALSE
		warning("Option summary.only=TRUE is not supported by Pbugs. summary.only has been coerced to FALSE\n")
	}
 
	if (is.R()) .fileCopy <- file.copy
	else .fileCopy <- splus.file.copy

	#Creates, and copies WinBUGS copies to, pbugs.directory 
	if(file.exists(pbugs.directory)){
		for(i in 1:n.chains){
			exists<-file.exists(paste(file.path(pbugs.directory,basename(bugs.directory)),"-",i,sep=""))
			if(!exists){
				isok1<-.fileCopy(bugs.directory,pbugs.directory,recursive=TRUE)
				isok2<-file.rename(file.path(pbugs.directory,basename(bugs.directory)),file.path(pbugs.directory,paste(basename(bugs.directory),"-",i,sep="")))
				if(!isok1 | !isok2){stop("Cannot create WinBUGS/OpenBUGS copies")}
			}
		}
	}else{
		isok<-dir.create(pbugs.directory,recursive=TRUE)
		if(!isok){stop(paste("Cannot create directory:",pbugs.directory,"\n"))}
		for(i in 1:n.chains){
			isok1<-.fileCopy(bugs.directory,pbugs.directory,recursive=TRUE)
			isok2<-file.rename(file.path(pbugs.directory,basename(bugs.directory)),file.path(pbugs.directory,paste(basename(bugs.directory),"-",i,sep="")))
			if(!isok1 | !isok2){stop("Cannot create OpenBUGS copies")}
		}
	}
	#
	#####################
	
	if (!is.null(working.directory)) {
		working.directory <- path.expand(working.directory)
		savedWD <- getwd()
		setwd(working.directory)
		on.exit(setwd(savedWD))
	}
	program <- match.arg(program)
	if (missing(bugs.directory) && !is.null(bugs.dir <- getOption("R2WinBUGS.bugs.directory"))) {
		bugs.directory <- bugs.dir
	}
	#if (program %in% c("openbugs", "OpenBUGS", "OpenBugs")) {
	#	if (!is.R()) stop("OpenBUGS is not yet available in S-PLUS")
	#	return(openbugs(data, inits, parameters.to.save, model.file, n.chains, n.iter, n.burnin, n.thin, n.sims, DIC = DIC, bugs.directory, working.directory, digits, over.relax = over.relax, seed = bugs.seed))
	#}
	if (!missing(inits) && !is.function(inits) && !is.null(inits) && (length(inits) != n.chains)) stop("Number of initialized chains (length(inits)) != n.chains")
	if (useWINE) {
		if (!is.R()) stop("Non-Windows platforms not yet supported in R2WinBUGS for S-PLUS")
		if (is.null(WINE)) WINE <- findUnixBinary(x = "wine")
		if (is.null(WINEPATH)) WINEPATH <- findUnixBinary(x = "winepath")
	}
	inTempDir <- FALSE
	if (is.null(working.directory)) {
		working.directory <- tempdir()
		if (useWINE) {
			working.directory <- gsub("//", "/", working.directory)
			Sys.chmod(working.directory, mode = "770")
			on.exit(Sys.chmod(working.directory, mode = "700"), add = TRUE)
		}
		savedWD <- getwd()
		setwd(working.directory)
		on.exit(setwd(savedWD), add = TRUE)
		inTempDir <- TRUE
	}
	if (is.function(model.file)) {
		temp <- tempfile("model")
		temp <- if (is.R() || .Platform$OS.type != "windows") {paste(temp, "txt", sep = ".")} else {gsub("\\.tmp$", ".txt", temp)}
		write.model(model.file, con = temp, digits = digits)
		model.file <- gsub("\\\\", "/", temp)
		if (!is.R()) on.exit(file.remove(model.file), add = TRUE)
	}
	#if (inTempDir && basename(model.file) == model.file) try(file.copy(file.path(savedWD, model.file), model.file, overwrite = TRUE))
 
	#####################
	# Pbugs-specific code
	if (inTempDir && basename(model.file) == model.file) try(.fileCopy(file.path(savedWD, model.file), model.file, overwrite = TRUE))
	#
	#####################
 
	if (!file.exists(model.file)) stop(paste(model.file, "does not exist."))
	if (file.info(model.file)$isdir) stop(paste(model.file, "is a directory, but a file is required."))
	if (!(length(data) == 1 && is.vector(data) && is.character(data) && (regexpr("\\.txt$", data) > 0))) {
		bugs.data.file <- bugs.data(data, dir = getwd(), digits)
	}
	else {
		#if (inTempDir && all(basename(data) == data)) try(file.copy(file.path(savedWD, data), data, overwrite = TRUE))
 
		#####################
		# Pbugs-specific code
		if (inTempDir && all(basename(data) == data)) try(.fileCopy(file.path(savedWD, data), data, overwrite = TRUE))
		#
		#####################
 
		if (!file.exists(data)) stop("File", data, "does not exist.")
		bugs.data.file <- data
	}
	
	if (is.character(inits)) {
		#if (inTempDir && all(basename(inits) == inits)) try(file.copy(file.path(savedWD, inits), inits, overwrite = TRUE))
 
		#####################
		# Pbugs-specific code
		if (inTempDir && all(basename(inits) == inits)) try(.fileCopy(file.path(savedWD, inits), inits, overwrite = TRUE))
		#
		#####################

		if (!all(file.exists(inits))) stop("One or more inits files are missing")
		if (length(inits) != n.chains) stop("Need one inits file for each chain")
		bugs.inits.files <- inits
	}
	else {
		if (!is.function(inits) && !is.null(inits) && (length(inits) != n.chains)) stop("Number of initialized chains (length(inits)) != n.chains")
		bugs.inits.files <- R2WinBUGS:::bugs.inits(inits, n.chains, digits)
	}
	if (DIC) parameters.to.save <- c(parameters.to.save, "deviance")
	if (!length(grep("\\.txt$", tolower(model.file)))) {
		new.model.file <- paste(basename(model.file), ".txt", sep = "")
		if (!is.null(working.directory)) new.model.file <- file.path(working.directory, new.model.file)
		#file.copy(model.file, new.model.file, overwrite = TRUE)

		#####################
		# Pbugs-specific code
		.fileCopy(model.file, new.model.file, overwrite = TRUE)
		#
		#####################
 
		on.exit(try(file.remove(new.model.file)), add = TRUE)
	}
	else {new.model.file <- model.file}
	if (useWINE) {
		new.model.file <- gsub("//", "/", new.model.file)
	}
	 
	#####################
	# Pbugs-specific code
	try(dir.create(file.path(working.directory,"Pbugs-working"),showWarnings=F))
	for(i in 1:n.chains){
		working.aux<-file.path(working.directory,"Pbugs-working",paste("ch",i,sep=""))
		try(dir.create(working.aux,showWarnings=F))
		try(.fileCopy(basename(model.file),file.path(working.aux,basename(model.file)),overwrite=TRUE))
		try(.fileCopy("data.txt",file.path(working.aux,"data.txt"),overwrite=TRUE))
		try(.fileCopy(paste("inits",i,".txt",sep=""),file.path(working.aux,"inits1.txt"),overwrite=TRUE))
 		setwd(working.aux)
  	R2WinBUGS:::bugs.script(parameters.to.save, 1, n.iter, n.burnin, n.thin, new.model.file, debug = debug, is.inits = !is.null(inits), bin = bin, DIC = DIC, useWINE = useWINE, newWINE = newWINE, WINEPATH = WINEPATH, bugs.seed = bugs.seed, summary.only = summary.only, save.history = save.history, bugs.data.file = bugs.data.file, bugs.inits.files = bugs.inits.files[1], over.relax = over.relax)
#		aux<-readLines("script.txt")
##		aux2<-sapply(aux,function(x){gsub(pattern=working.directory,x=x,replacement=file.path(working.directory,"Pbugs",paste("ch",i,sep="")))})
##		writeLines(aux,file.path(working.directory,"Pbugs",paste("ch",i,sep=""),"script.txt"))
		setwd(working.directory)
	}
 
	PWinBUGS.run(n.burnin, bugs.directory, cluster, pbugs.directory, n.chains, WINE = WINE, useWINE = useWINE, newWINE = newWINE, WINEPATH = WINEPATH)
 
	#
	#####################
 
	if (codaPkg) return(file.path(getwd(), paste("coda", 1:n.chains, ".txt", sep = "")))
	#if (summary.only) return(bugs.log("log.txt"))
	#sims <- c(R2WinBUGS:::bugs.sims(parameters.to.save, n.chains, n.iter, n.burnin, n.thin, DIC=DIC), model.file = model.file, program = program)
	#####################
	# Pbugs-specific code (modified from R2WinBUGS:::bugs.sims)
	sims <- c(R2WinBUGS:::bugs.sims(parameters.to.save, n.chains, n.iter, n.burnin, n.thin, DIC=FALSE), model.file = model.file, program = program)
	if (DIC) {
		LOG<-list(length=n.chains)
		pD<-rep(NA,n.chains)
		DIC<-rep(NA,n.chains)
		for(i in 1:n.chains) LOG[[i]] <- R2WinBUGS:::bugs.log(file.path(working.directory,"Pbugs-working",paste("ch",i,sep=""),"log.txt"))$DIC
		if (any(is.na(LOG))) {
			deviance <- sims$sims.array[, , dim(sims.array)[3], drop = FALSE]
			if (!is.R()) dimnames(deviance) <- NULL
			dim(deviance) <- dim(deviance)[1:2]
			pD <- numeric(n.chains)
			DIC <- numeric(n.chains)
			for (i in 1:n.chains) {
				pD[i] <- var(deviance[, i])/2
				DIC[i] <- mean(deviance[, i]) + pD[i]
			}
			sims$DICbyR<-TRUE
		}
		else {
			for(i in 1:n.chains){
				pD[i]<-LOG[[i]][nrow(LOG[[i]]), 3]
				DIC[i]<-LOG[[i]][nrow(LOG[[i]]), 4]
			}
			sims$DICbyR<-FALSE
		}
		sims$isDIC<-TRUE
		sims$pD<-mean(pD)
		sims$DIC<-mean(DIC)
	}
	#
	#####################
 
	#if (clearWD) file.remove(c(bugs.data.file, "log.odc", "log.txt", "codaIndex.txt", bugs.inits.files, "script.txt", paste("coda", 1:n.chains, ".txt", sep = "")))
 
	#####################
	# Pbugs-specific code
	if (clearWD) {
		### COMPROBAR QUE FUNCIONA BIEN ###
		file.remove(c(bugs.data.file, "log.odc", "log.txt", "codaIndex.txt", bugs.inits.files, "script.txt", paste("coda", 1:n.chains, ".txt", sep = ""),"Pbugs"))
	}
	#
	#####################
 
	#class(sims) <- c("WinBUGS","bugs")
	class(sims) <- "bugs"
	return(sims)
}
 
 
 
PWinBUGS.run<-function (n.burnin, bugs.directory, cluster, pbugs.directory, n.chains, useWINE = .Platform$OS.type != "windows", WINE = NULL, newWINE = TRUE, WINEPATH = NULL){
	if (useWINE && !is.R()) stop("Non-Windows platforms not yet supported in R2WinBUGS for S-PLUS")
	if (useWINE && (substr(bugs.directory, 2, 2) == ":")) bugs.directory <- R2WinBUGS:::win2native(bugs.directory, newWINE = newWINE, WINEPATH = WINEPATH)
	#try(R2WinBUGS:::bugs.update.settings(n.burnin, bugs.directory))
	#####################
	# Pbugs-specific code
	for(i in 1:n.chains){
	  try(R2WinBUGS:::bugs.update.settings(n.burnin, file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""))))
	}
	#####################
	if (is.R()) .fileCopy <- file.copy
	else .fileCopy <- splus.file.copy
	#on.exit(IsItOk<-try(.fileCopy(file.path(pbugs.directory,"System/Rsrc/Registry_Rsave.odc"), file.path(pbugs.directory, "System/Rsrc/Registry.odc"), overwrite = TRUE)))
	#on.exit(while(!IsItOk & counter<6){IsItOk<-try(.fileCopy(file.path(pbugs.directory, "System/Rsrc/Registry_Rsave.odc"), file.path(pbugs.directory, "System/Rsrc/Registry.odc"), overwrite = TRUE));Sys.sleep(0.5);counter<-counter+1},add=TRUE)
	#####################
	# Pbugs-specific code
	#for(i in 1:n.chains){
	#  on.exit(try(.fileCopy(file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""),"System/Rsrc/Registry_Rsave.odc"), file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""), "System/Rsrc/Registry.odc"), overwrite = TRUE)), add = TRUE)
	  #on.exit(IsItOk<-try(.fileCopy(file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""),"System/Rsrc/Registry_Rsave.odc"), file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""), "System/Rsrc/Registry.odc"), overwrite = TRUE)))
	  #counter<-1
	  #on.exit(while(!IsItOk & counter<6){IsItOk<-try(.fileCopy(file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""),"System/Rsrc/Registry_Rsave.odc"), file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""),"System/Rsrc/Registry.odc"), overwrite = TRUE));Sys.sleep(0.5);counter<-counter+1},add=TRUE)
	#}
	on.exit(try(.fileCopy(file.path(pbugs.directory, paste(basename(bugs.directory),"-",1:3,sep=""),"System/Rsrc/Registry_Rsave.odc"), file.path(pbugs.directory, paste(basename(bugs.directory),"-",1:3,sep=""), "System/Rsrc/Registry.odc"), overwrite = TRUE)), add = TRUE)
	  #on.exit(IsItOk<-try(.fileCopy(file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""),"System/Rsrc/Registry_Rsave.odc"), file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""), "System/Rsrc/Registry.odc"), overwrite = TRUE)))
	  #counter<-1
	  #on.exit(while(!IsItOk & counter<6){IsItOk<-try(.fileCopy(file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""),"System/Rsrc/Registry_Rsave.odc"), file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""),"System/Rsrc/Registry.odc"), overwrite = TRUE));Sys.sleep(0.5);counter<-counter+1},add=TRUE)
	
	#
	#####################
	#dos.location <- file.path(bugs.directory, grep("^Win[[:alnum:]]*[.]exe$", list.files(bugs.directory), value = TRUE)[1])
	#if (!file.exists(dos.location)) stop(paste("WinBUGS executable does not exist in", bugs.directory))

	#####################
	# Pbugs-specific code
	bugsCall<-vector(length=n.chains)
	for(i in 1:n.chains){
		dos.location <- file.path(pbugs.directory, paste(basename(bugs.directory),"-",i,sep=""),grep("^Win[[:alnum:]]*[.]exe$", list.files(bugs.directory), value = TRUE)[1])
		if (!file.exists(dos.location)) stop(paste("WinBUGS executable does not exist in", paste(basename(bugs.directory),"-",i,sep="")))
		bugsCall[i] <- paste("\"", dos.location, "\" /par \"", R2WinBUGS:::native2win(file.path(getwd(),"Pbugs-working", paste("ch",i,sep=""), "script.txt"), useWINE = useWINE, newWINE = newWINE, WINEPATH = WINEPATH), "\"", sep = "")
    		if (useWINE) bugsCall[i] <- paste(WINE, bugsCall[i])
		#cat(i, bugsCall[i],"\n")
	}
	if(is.null(cluster)){
		cl<<-makeCluster(n.chains,type="SOCK")
		on.exit(stopCluster(cl),add=TRUE)
	}else{
		cl<<-cluster
	}
	#cat("Virtual cluster initiated","\n")
	clusterEvalQ(cl,library(R2WinBUGS))
	clusterEvalQ(cl,set.seed(1))
	temp<-clusterApply(cl,bugsCall,system)
	.fileCopy(file.path(getwd(), "Pbugs-working", "ch1", "codaIndex.txt"),"codaIndex.txt",overwrite=TRUE)
	for(i in 1:n.chains) .fileCopy(file.path(getwd(), "Pbugs-working", paste("ch",i,sep=""), "coda1.txt"),paste("coda",i,".txt",sep=""),overwrite=TRUE)
	#if (temp == -1) stop("Error in bugs.run().\nCheck that WinBUGS is in the specified directory.")
	if (any(unlist(temp) == -1)) stop("Error in bugs.run().\nCheck that WinBUGS is in the specified directory.")
	#
	#####################
 
	if (is.R()) tmp <- scan("coda1.txt", character(), quiet = TRUE, sep = "\n")
	else tmp <- scan("coda1.txt", character(), sep = "\n")
	if (length(grep("BUGS did not run correctly", tmp)) > 0) stop("Look at the log file and\ntry again with 'debug=TRUE' to figure out what went wrong within Bugs.")
}
 
 
 
POpenBUGS<-function (data, inits, parameters.to.save, n.iter, model.file = "model.txt", n.chains = 3, n.burnin = floor(n.iter/2), n.thin = 1, saveExec = FALSE, restart = FALSE, debug = FALSE, DIC = TRUE, digits = 5, codaPkg = FALSE, OpenBUGS.pgm = NULL, working.directory = NULL, clearWD = FALSE, useWINE = FALSE, WINE = NULL, newWINE = TRUE, WINEPATH = NULL, bugs.seed = 1, summary.only = FALSE, save.history = (.Platform$OS.type == "windows" | useWINE == TRUE), over.relax = FALSE,cluster, pbugs.directory) {
 
	#####################
	# Pbugs-specific code
	if(summary.only) {
		summary.only<-FALSE
		warning("Option summary.only=TRUE is not supported by Pbugs. summary.only has been coerced to FALSE\n")
	}
  
  if (is.null(OpenBUGS.pgm)) {
    OpenBUGS.pgm <- R2OpenBUGS:::findOpenBUGS()
    if (.Platform$OS.type == "windows" || useWINE) OpenBUGS.pgm <- file.path(OpenBUGS.pgm, "OpenBUGS.exe")
  }
  else if (OpenBUGS.pgm == "OpenBUGS") OpenBUGS.pgm <- Sys.which("OpenBUGS")
  if (!file.exists(OpenBUGS.pgm)) stop("Cannot find the OpenBUGS program")
  if (useWINE && (substr(OpenBUGS.pgm, 2, 2) == ":")) OpenBUGS.pgm <- win2native(OpenBUGS.pgm, newWINE = newWINE, WINEPATH = WINEPATH)

  #####################
  # Pbugs-specific code
	# Creates, and copies WinBUGS copies to, pbugs.directory 
	OpenBUGS.path<-dirname(OpenBUGS.pgm)
	OpenBUGS.folder<-strsplit(OpenBUGS.path,"/")[[1]][length(strsplit(OpenBUGS.path,"/")[[1]])]
	if(file.exists(pbugs.directory)){
		for(i in 1:n.chains){
			exists<-file.exists(paste(file.path(pbugs.directory,OpenBUGS.folder),"-",i,sep=""))
			if(!exists){
				isok1<-file.copy(OpenBUGS.path,pbugs.directory,recursive=TRUE)
				isok2<-file.rename(file.path(pbugs.directory,OpenBUGS.folder),file.path(pbugs.directory,paste(OpenBUGS.folder,"-",i,sep="")))
				if(!isok1 | !isok2){stop("Cannot create OpenBUGS copies")}
			}
		}
	}else{
	  isok<-dir.create(pbugs.directory,recursive=TRUE)
		if(!isok){stop(paste("Cannot create directory: ",pbugs.directory,"\n"))}
		for(i in 1:n.chains){
			isok1<-file.copy(OpenBUGS.path,pbugs.directory,recursive=TRUE)
			isok2<-file.rename(file.path(pbugs.directory,OpenBUGS.folder),file.path(pbugs.directory,paste(OpenBUGS.folder,"-",i,sep="")))
			if(!isok1 | !isok2){stop("Cannot create OpenBUGS copies")}
		}
	}
	#
	#####################
  if (.Platform$OS.type != "windows" && !useWINE) {
		if (debug) stop("The debug option is not available with linux/unix")
		if (save.history) stop("History plots (save.history) are not available with linux/unix")
	}
	if (!bugs.seed %in% 1:14) stop("OpenBUGS seed must be integer in 1:14")
	if (!is.function(model.file) && length(grep("\\.bug", tolower(model.file)))) stop("model.file must be renamed with .txt rather than .bug")
	if (is.null(working.directory) && (saveExec || restart)) stop("The working directory must be specified when saveExec or restart is TRUE")
	if (!is.null(working.directory)) {
		working.directory <- path.expand(working.directory)
		savedWD <- getwd()
		setwd(working.directory)
		on.exit(setwd(savedWD))
	}
	if (!missing(inits) && !is.function(inits) && !is.null(inits) && (length(inits) != n.chains)) stop("Number of initialized chains (length(inits)) != n.chains")
	if (useWINE) {
		if (is.null(WINE)) WINE <- findUnixBinary(x = "wine")
		if (is.null(WINEPATH)) WINEPATH <- findUnixBinary(x = "winepath")
	}
	inTempDir <- FALSE
	if (is.null(working.directory)) {
		working.directory <- tempdir()
		if (useWINE) {
			working.directory <- gsub("//", "/", working.directory)
			Sys.chmod(working.directory, mode = "770")
			on.exit(Sys.chmod(working.directory, mode = "700"), add = TRUE)
		}
		savedWD <- getwd()
		setwd(working.directory)
		on.exit(setwd(savedWD), add = TRUE)
		inTempDir <- TRUE
	}
	if (is.function(model.file)) {
		temp <- tempfile("model")
		temp <- paste(temp, "txt", sep = ".")
		write.model(model.file, con = temp, digits = digits)
		model.file <- gsub("\\\\", "/", temp)
	}
	if (inTempDir && basename(model.file) == model.file) try(file.copy(file.path(savedWD, model.file), model.file, overwrite = TRUE))
	if (!file.exists(model.file)) stop(paste(model.file, "does not exist."))
	if (file.info(model.file)$isdir) stop(paste(model.file, "is a directory, but a file is required."))
	if (!(length(data) == 1 && is.vector(data) && is.character(data) && (regexpr("\\.txt$", data) > 0))) bugs.data.file <- bugs.data(data, dir = getwd(), digits)
	else {
		if (inTempDir && all(basename(data) == data)) try(file.copy(file.path(savedWD, data), data, overwrite = TRUE))
		if (!file.exists(data)) stop("File", data, "does not exist.")
		bugs.data.file <- data
	}
	if (is.character(inits)) {
		if (inTempDir && all(basename(inits) == inits)) try(file.copy(file.path(savedWD, inits), inits, overwrite = TRUE))
		if (!all(file.exists(inits))) stop("One or more inits files are missing")
    if (length(inits) != n.chains) stop("Need one inits file for each chain")
		bugs.inits.files <- inits
	}
	else {
		if (!is.function(inits) && !is.null(inits) && (length(inits) != n.chains)) stop("Number of initialized chains (length(inits)) != n.chains")
		bugs.inits.files <- R2OpenBUGS:::bugs.inits(inits, n.chains, digits)
	}
	if (DIC) parameters.to.save <- c(parameters.to.save, "deviance")
	if (!length(grep("\\.txt$", tolower(model.file)))) {
		new.model.file <- paste(basename(model.file), ".txt", sep = "")
		if (!is.null(working.directory)) new.model.file <- file.path(working.directory, new.model.file)
		file.copy(model.file, new.model.file, overwrite = TRUE)
		on.exit(try(file.remove(new.model.file)), add = TRUE)
	}
	else new.model.file <- model.file
	model.file.bug <- gsub("\\.txt", ".bug", basename(new.model.file))
	if (restart && !file.exists(model.file.bug)) stop("The .bug restart file was not found in the working directory")
	if (useWINE) new.model.file <- gsub("//", "/", new.model.file)

	#####################
	# Pbugs-specific code
	try(dir.create(file.path(working.directory,"Pbugs-working"),showWarnings=F))
	for(i in 1:n.chains){
	  working.aux<-file.path(working.directory,"Pbugs-working",paste("ch",i,sep=""))
		try(dir.create(working.aux,showWarnings=F))
		try(file.copy(basename(model.file),file.path(working.aux,basename(model.file)),overwrite=TRUE))
		try(file.copy("data.txt",file.path(working.aux,"data.txt"),overwrite=TRUE))
		try(file.copy(paste("inits",i,".txt",sep=""),file.path(working.aux,"inits1.txt"),overwrite=TRUE))
	  setwd(working.aux)
	  R2OpenBUGS:::bugs.script(parameters.to.save, 1, n.iter, n.burnin, n.thin, saveExec, restart, model.file.bug, new.model.file, debug = debug, is.inits = !is.null(inits), DIC = DIC, useWINE = useWINE, newWINE = newWINE, WINEPATH = WINEPATH, bugs.seed = bugs.seed, summary.only = summary.only, save.history = save.history, bugs.data.file = bugs.data.file, bugs.inits.files = bugs.inits.files[1], over.relax = over.relax)
	  setwd(working.directory)
	}
	#R2OpenBUGS:::bugs.script(parameters.to.save, 1, n.iter, n.burnin, n.thin, saveExec, restart, model.file.bug, new.model.file, debug = debug, is.inits = !is.null(inits), DIC = DIC, useWINE = useWINE, newWINE = newWINE, WINEPATH = WINEPATH, bugs.seed = bugs.seed, summary.only = summary.only, save.history = save.history, bugs.data.file = bugs.data.file, bugs.inits.files = bugs.inits.files[1], over.relax = over.relax)
	#aux<-readLines("script.txt")
	#for(i in 1:n.chains){
	#	aux2<-sapply(aux,function(x){gsub(pattern=working.directory,x=x,replacement=file.path(working.directory,"Pbugs",paste("ch",i,sep="")))})
	#	writeLines(aux2,file.path(working.directory,"Pbugs",paste("ch",i,sep=""),"script.txt"))
	#}
 
  #bugs.run(n.burnin, OpenBUGS.pgm, debug = debug, WINE = WINE, useWINE = useWINE, newWINE = newWINE, WINEPATH = WINEPATH)
	POpenBUGS.run(n.burnin, OpenBUGS.pgm, cluster, pbugs.directory, n.chains, debug=debug, WINE = WINE, useWINE = useWINE, newWINE = newWINE, WINEPATH = WINEPATH)
	
	#
	#####################
 
	if (codaPkg) return(file.path(getwd(), paste("CODAchain", 1:n.chains, ".txt", sep = "")))
	#if (summary.only) return(bugs.log("log.txt"))
	#sims <- c(R2OpenBUGS:::bugs.sims(parameters.to.save, n.chains, n.iter, n.burnin, n.thin, DIC=FALSE), model.file = model.file)
	#	if (clearWD) file.remove(c(bugs.data.file, "log.odc", "log.txt", "CODAIndex.txt", bugs.inits.files, "script.txt", paste("CODAchain", 1:n.chains, ".txt", sep = "")))
 
	#####################
	# Pbugs-specific code (modified from R2WinBUGS:::bugs.sims)
  sims <- c(R2OpenBUGS:::bugs.sims(parameters.to.save, n.chains, n.iter, n.burnin, n.thin, DIC=FALSE), model.file = model.file)
	if (DIC) {
		LOG<-list(length=n.chains)
		pD<-rep(NA,n.chains)
		DIC<-rep(NA,n.chains)
		for(i in 1:n.chains) LOG[[i]] <- R2OpenBUGS:::bugs.log(file.path(working.directory,"Pbugs-working",paste("ch",i,sep=""),"log.txt"))$DIC
		if (any(is.na(LOG))) {
			deviance <- sims$sims.array[, , dim(sims.array)[3], drop = FALSE]
			if (!is.R()) dimnames(deviance) <- NULL
			dim(deviance) <- dim(deviance)[1:2]
			pD <- numeric(n.chains)
			DIC <- numeric(n.chains)
			for (i in 1:n.chains) {
				pD[i] <- var(deviance[, i])/2
				DIC[i] <- mean(deviance[, i]) + pD[i]
			}
			sims$DICbyR<-TRUE
		}
		else {
			for(i in 1:n.chains){
				pD[i]<-LOG[[i]][nrow(LOG[[i]]), 4]
				DIC[i]<-LOG[[i]][nrow(LOG[[i]]), 3]
			}
			sims$DICbyR<-FALSE
		}
		sims$isDIC<-TRUE
		sims$pD<-mean(pD)
		sims$DIC<-mean(DIC)
	}
	#
	#####################
 
	#####################
	# Pbugs-specific code
	if (clearWD) {
		#file.remove(c(bugs.data.file, "log.odc", "log.txt", "codaIndex.txt", bugs.inits.files, "script.txt", paste("coda", 1:n.chains, ".txt", sep = "")))
		### COMPROBAR QUE FUNCIONA BIEN ###
		file.remove(c(bugs.data.file, "log.odc", "log.txt", "codaIndex.txt", bugs.inits.files, "script.txt", paste("coda", 1:n.chains, ".txt", sep = "")),"Pbugs-working")
	}
  #class(sims) <- c("OpenBUGS","bugs")
  #
	#####################
 
  class(sims) <- "bugs"
	sims
}
 
 
POpenBUGS.run<-function (n.burnin, OpenBUGS.pgm, cluster, pbugs.directory, n.chains, debug = FALSE, useWINE = FALSE, WINE = NULL, newWINE = TRUE, WINEPATH = NULL){
	#if (.Platform$OS.type == "windows" || useWINE) {
	#	bugsCall <- paste("\"", OpenBUGS.pgm, "\" /PAR \"", native2win(file.path(getwd(), "script.txt"), useWINE = useWINE, newWINE = newWINE, WINEPATH = WINEPATH), "\" /", sep = "")
	#	if (!debug) bugsCall <- paste(bugsCall, "HEADLESS", sep = "")
	#	if (useWINE) bugsCall <- paste(WINE, bugsCall)
	#}
	#else bugsCall <- paste(OpenBUGS.pgm, "<", "script.txt", ">", file.path(getwd(), "log.txt"))
 
	#####################
	# Pbugs-specific code	
	OpenBUGS.path<-dirname(OpenBUGS.pgm)
	OpenBUGS.folder<-strsplit(OpenBUGS.path,"/")[[1]][length(strsplit(OpenBUGS.path,"/")[[1]])]

	bugsCall<-vector(length=n.chains)
	if (.Platform$OS.type == "windows" || useWINE) {
		for(i in 1:n.chains){
		  Openfile<-file.path(paste(file.path(pbugs.directory,OpenBUGS.folder),"-",i,sep=""),basename(OpenBUGS.pgm))
			bugsCall[i] <- paste("\"", Openfile, "\" /PAR \"", R2OpenBUGS:::native2win(file.path(getwd(), "Pbugs-working", paste("ch",i,sep=""), "script.txt"), useWINE = useWINE, newWINE = newWINE, WINEPATH = WINEPATH), "\" /", sep = "")
			#Next sentence avoids openbugs to be launched (originally en R2OpenBUGS).
			#if (!debug) bugsCall[i] <- paste(bugsCall[i], "HEADLESS", sep = "")
			if (useWINE) bugsCall[i] <- paste(WINE, bugsCall[i])
		}
	}
	else{ 
	  for(i in 1:n.chains){ 
	    Openfile<-file.path(paste(file.path(pbugs.directory,OpenBUGS.folder),"-",i,sep=""),basename(OpenBUGS.pgm))
	    bugsCall[i] <- paste(Openfile, "<", file.path(getwd(),"Pbugs-working", paste("ch",i,sep=""),"script.txt"), ">", file.path(getwd(), "Pbugs-working", paste("ch",i,sep=""), "log.txt"))
	  }
	}
	if(is.null(cluster)){
	  cl<<-makeCluster(n.chains,type="SOCK")
	  on.exit(stopCluster(cl),add=TRUE)
	}else{
	  cl<<-cluster
	}
	#cat("Virtual cluster initiated","\n")
	clusterEvalQ(cl,library(R2OpenBUGS))
	if ((.Platform$OS.type == "windows" || useWINE) && debug) temp<-clusterApply(cl,bugsCall,system, invisible=FALSE)
	else temp<-clusterApply(cl,bugsCall,system)
	file.copy(file.path(getwd(), "Pbugs-working", "ch1", "CODAindex.txt"),"CODAindex.txt",overwrite=TRUE)
	for(i in 1:n.chains) file.copy(file.path(getwd(), "Pbugs-working", paste("ch",i,sep=""), "CODAchain1.txt"),paste("CODAchain",i,".txt",sep=""),overwrite=TRUE)
	if (any(unlist(temp) == -1)) stop("Error in bugs.run(). Check that OpenBUGS is in the specified directory.")
	tmp<-list(length=n.chains)
	for(i in 1:n.chains){
		tmp[[i]] <- scan(paste("CODAchain",i,".txt",sep=""), character(), quiet = TRUE, sep = "\n")
		tmp[[i]] <- tmp[[i]][1:min(100, length(tmp[[i]]))]
		tmp[[i]]<-grep("OpenBUGS did not run correctly", tmp[[i]])
	}
	if (any(sapply(tmp,length) > 0)) stop(paste("Look at the log file in ", getwd(), " and\ntry again with 'debug=TRUE' to figure out what went wrong within OpenBUGS."))
	#
	#####################
}

#findUnixBinary <- function(x)
#{
#  ## --- Environmental variable ---
#  tmp <- Sys.getenv(toupper(x))
#  if(nchar(tmp) != 0 && file.exists(tmp)) return(tmp)
#  ## else
#
#  ## --- Standard place ---
#  tmp <- paste("/usr/bin", x, sep="")
#  if(file.exists(tmp)) return(tmp)
#  ## else ...
#
#  ## --- Which ---
#  tmp <- system(paste("which ", x, sep=""), intern=TRUE)
#  if(length(tmp) != 0 && file.exists(tmp)) return(tmp)
#  ## else ..
#
#  ## --- Locate ---
#  tmp <- system(paste("locate ", x, " | grep bin/", x, "$", sep=""), intern=TRUE)
#  tmp <- tmp[length(tmp)] ## keep only last hit
#  if(length(tmp) > 0 && file.exists(tmp)) return(tmp)
#
#  stop(paste("couldn't find", x, "binary file"))
#}

