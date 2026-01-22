compute_dynamic_value <- function(landings, min_price_age, max_price_age, breakpoints=c(15, 30)){
   if(landings < breakpoints[1]){
      return(max_price_age)
   }else if(landings >= breakpoints[1] & landings <= breakpoints[2]){
      return(
            min_price_age + (breakpoints[2]-landings)/(breakpoints[2]-breakpoints[1])*(max_price_age-min_price_age)
      )
   }else{
      return(min_price_age)
   }
}


price_age_f_low <- c(0.597895623, 1.320303448, 1.320303448, 1.856562267, 2.610111345, 2.610111345, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875)
price_age_m_low <- c(0.597895623, 0.597895623, 1.320303448, 1.320303448, 1.856562267, 1.856562267, 1.856562267, 1.856562267, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345)
price_data_low <- matrix(c(price_age_f_low, price_age_m_low), nrow=length(price_age_f_low), ncol=2)
dimnames(price_data_low) <- list("age"=2:31, "sex"=c("F", "M"))

price_age_f_max <- c(7.917460094, 8.40756497, 8.40756497, 9.944657109, 11.46480347, 11.46480347, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658)
price_age_m_max <- c(7.917460094, 7.917460094, 8.40756497, 8.40756497, 9.944657109, 9.944657109, 9.944657109, 9.944657109, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347)
price_data_max <- matrix(c(price_age_f_max, price_age_m_max), nrow=length(price_age_f_max), ncol=2)
dimnames(price_data_max) <- list("age"=2:31, "sex"=c("F", "M"))


price_data_low %>% as_tibble() %>%
   mutate(age=2:31)


compute_dynamic_value(20, price_age_f_low, price_age_f_max)

catches <- seq(0, 40, by=1)

price_matrix_f <- sapply(catches, function(x) compute_dynamic_value(x, price_age_f_low, price_age_f_max))
price_matrix_m <- sapply(catches, function(x) compute_dynamic_value(x, price_age_m_low, price_age_m_max))

price_df_f <- as_tibble(price_matrix_f) %>%
   mutate(age=2:31) %>%
   pivot_longer(cols=-age, names_to="catch", values_to="price") %>%
   mutate(catch=as.integer(gsub("V", "", catch))-1) %>%
   mutate(
      sex = "F",
      price_grade = case_when(
         age %in% c(2) ~ "1 (Age 2)",
         age %in% c(3,4) ~ "2 (Age 3-4)",
         age %in% c(5) ~ "3 (Age 5)",
         age %in% c(6,7) ~ "4 (Age 6-7)",
         age %in% seq(8,14,1) ~ "5 (Age 8-14)",
         age %in% seq(15,31,1) ~ "6 (Age 15-31)"
      )
   )

price_df_m <- as_tibble(price_matrix_m) %>%
   mutate(age=2:31) %>%
   pivot_longer(cols=-age, names_to="catch", values_to="price") %>%
   mutate(catch=as.integer(gsub("V", "", catch))-1) %>%
   mutate(
      sex = "F",
      price_grade = case_when(
         age %in% c(2,3) ~ "1 (Age 2)",
         age %in% c(4,5) ~ "2 (Age 3-4)",
         age %in% c(6,7,8,9) ~ "3 (Age 5)",
         age %in% c(6,7) ~ "4 (Age 6-7)"
      )
   )


ggplot(price_df, aes(x=catch, y=price, color=factor(price_grade))) +
   geom_line() +
   labs(x="Landings (1000s mt)", y="Price per kilogram", color="Price Grade") +
   custom_theme
ggsave(file.path(here::here(), "figures", "hybrid_hcrs", "dynamic_price_curve.png"), width=8, height=6)
