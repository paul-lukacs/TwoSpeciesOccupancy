sites <- nrow( EHf )
occ <-  ncol( EHf )

for( i in 1:sites ){
	if( runif(1) < 0.6 ){
		EHf[i,] <- rbinom( nocc, 1, 0.05 )
	}
}

for( i in 1:sites ){
	if( runif(1) < 0.3 ){
		EHmf[i,] <- rbinom( nocc, 1, 0.1 )
	}
}