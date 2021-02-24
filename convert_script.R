# script to convert LCMS raw imputs into the concentration outputs

# libraries
library(tidyverse)
library(xlsx)
library(readxl)

##### useful functions #######

# Constructing Quadratic Formula
result_quad <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
    x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
    result = paste(x_1, x_2, sep=";")
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    x = -b/(2*a)
  }
  else {"There are no real roots."} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
  b^2-4*a*c
}

###### Parameters #######
mol_to_weight <- read_excel("molecule_tables.xlsx", sheet = "molecule_tables",
                            col_types=c("text", "numeric", "numeric", "numeric",
                                        "numeric","numeric","numeric"), na="NA")

alt_names <- read_excel("molecule_tables.xlsx", sheet = "alt_name",
                            col_types=c("text", "numeric", "numeric", "numeric")
                        , na="NA")

name_file <- "test_data.xlsx"

sample_list <- read_excel("sample_list.xlsx") %>% pull(ID)

##### processing samples data ##########

# import samples
for (sample_ID in sample_list){
    sample_name <- as.character(sample_ID)
    sample <- read_excel(name_file, sheet= sample_name) %>% separate(Chromatogram, into=c("EIC", "molecule_weight", "All", "MS"), sep=" ", remove=FALSE) %>%
      select(Chromatogram,`RT [min]`, Area, molecule_weight)
    
    # correct molecule weights
    sample$molecule_weight <- trunc(as.numeric(sample$molecule_weight))
    sample <- sample %>% mutate(molecule_weight=ifelse(molecule_weight==1217,608,molecule_weight))
    
    # find molecule name from weight
    sample_t1 <- left_join(sample,mol_to_weight,by = "molecule_weight")
    
    # normalized by standard
    stand <- sample_t1 %>% filter(molecule_weight==909) 
    sample_t2 <- sample_t1 %>% mutate(norm=Area/stand$Area)
    
    # solve quadratic
    sample_t3 <- sample_t2 %>% filter(!molecule_weight==909) %>% rowwise() %>% 
      mutate(quad=ifelse(!is.na(quad_a), result_quad(quad_a, quad_b, (quad_c-norm)), NA)) %>%
      separate(quad, into=c("quad_res1", "quad_res2"), sep=";")
    
    # solve linear
    sample_t3 <- sample_t3 %>% mutate(lin_res=((norm-lin_b)/lin_a))
    
    # normalization to 180ml*10
    sample_t4 <- sample_t3 %>% mutate(quad_res1=((as.numeric(quad_res1)/180)*200)*10) %>% 
      mutate(quad_res2=((as.numeric(quad_res2)/180)*200)*10) %>%
      mutate(lin_res=((as.numeric(lin_res)/180)*200)*10)
    
    # reshape and print results
    res_table <- sample_t4 %>% select(-quad_a, -quad_b, -quad_c, -lin_a, -lin_b, -norm)
    alt_table <- inner_join(res_table,alt_names) %>% rowwise() %>%
      filter(RT_start <= `RT [min]` && `RT [min]` <= RT_end) %>% 
      select(-RT_start, -RT_end)
    
    res_table <- left_join(res_table, alt_table)
    
    res_table <- data.frame(res_table)
    name_out <- paste("res", sample_ID, sep="")
    
    write.xlsx(res_table, "results.xlsx", sheetName = name_out, append=TRUE, row.names = FALSE, showNA=FALSE)
    
}




