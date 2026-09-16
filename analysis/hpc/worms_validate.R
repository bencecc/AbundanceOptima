#### ---- Utility function to validate fish taxonomies; used in combination with ---- #####
#### ---- fish_validate
###########################################################################################

worms_validate <- function(name) {
	
	try(verif <- wm_records_taxamatch(name)[[1]], 
			silent = T)
	
	if(exists("verif")) {
		try(aphia_h <- wm_classification(id = verif$AphiaID[1]), silent = T)
		if(exists("aphia_h")) {
			if(any(c("Actinopterygii","Elasmobranchii","Osteichthyes","Chondrichthyes") %in% aphia_h$scientificname)) {
				tmp.name <- as.character(verif$scientificname[1])
			}
		}
	}
	
	if(!exists("tmp.name")) {
		tmp.name <- data.frame()
	}
	
	as_tibble(tmp.name)
	
}
