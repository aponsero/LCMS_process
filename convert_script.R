# script to convert LCMS raw imputs into the concentration outputs

# libraries
library("readxl")
library("tidyverse")

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
mol_to_weight <- read_csv("molecule_tables.csv") 

sample_list=c("sample11.csv",
              "sample37.csv",
              "sample42.csv",
              "sample60.csv",
              "sample69.csv")

##### processing samples data ##########

# import samples
for (sample_file in sample_list){
    sample <- read_csv(sample_file) %>% separate(Chromatogram, into=c("EIC", "molecule_weight", "All", "MS"), sep=" ", remove=FALSE) %>%
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
      mutate(quad=result_quad(quad_a, quad_b, (quad_c-norm))) %>%
      separate(quad, into=c("quad_res1", "quad_res2"), sep=";")
    
    # solve linear
    sample_t3 <- sample_t3 %>% mutate(lin_res=((norm-lin_b)/lin_a))
    
    # normalization to 180ml*10
    sample_t4 <- sample_t3 %>% mutate(quad_res1=((as.numeric(quad_res1)/180)*200)*10) %>% 
      mutate(quad_res2=((as.numeric(quad_res2)/180)*200)*10) %>%
      mutate(lin_res=((as.numeric(lin_res)/180)*200)*10)
    
    # reshape and print results
    quad_table <- sample_t4 %>% select(-quad_a, -quad_b, -quad_c, -lin_a, -lin_b)
    out_file <- paste("res", sample_file, sep="_")
    write_csv(quad_table, out_file)
}




