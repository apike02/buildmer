calcor <- function (x) {
	x <- exp(x)
	if (x < 1) paste0('1:',custround(1/x,neg=F,trunc=T)     )
	else       paste0(     custround(  x,neg=F,trunc=T),':1')
}

custround <- function (i,neg=T,trunc=F) {
	i = as.numeric(i) #the matrix changes to a character matrix because of calcor
	prec <- 3
	ir <- round(i,prec)
	while (i > 0 && ir == 0) {
		prec <- prec + 1
		ir <- round(i,prec+1) #3 -> 5 -> 6 -> 7 -> ...
	}
	ir <- as.character(ir)
	while (nchar(sub('.*\\.','',ir)) < prec) ir <- paste0(ir,'0')
	if (trunc) ir <- sub('^0+','',ir)
	if (neg & i >= 0) ir <- paste0('\\hphantom{-}',ir)
	ir
}

nohp <- function (x) sub('\\hphantom{-}','',x,fixed=T)

paperify <- function (x,aliases) {
	if (!x %in% names(aliases)) return(x)
	x <- names(aliases) == x
	x <- unlist(aliases[x]) #x <- aliases[[x]] gives 'attempt to select less than one element' error??
}

stars <- function (x)  return(if (as.numeric(x) < .001) '$**$$*$'
			else if (as.numeric(x) < .01) '$**$'
			else if (as.numeric(x) < .05) '$*$'
			else '')

tblprintln <- function (x) {
	l <- paste0(x,collapse=' & ')
	cat(l,'\\\\\n',sep='')
}
