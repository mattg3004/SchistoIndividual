data_2000 = read_excel("data_gurarie_milalani.xlsx", sheet = "Age and egg count Milalani 2000")
data_2000$Age = data_2000$AGE
data_2000 = data_2000[,-2]

data_2003 = read_excel("data_gurarie_milalani.xlsx", sheet = "Milalani 2003")

data_2009 = read_excel("data_gurarie_milalani.xlsx", sheet = "Age_and_Egg_count_Milalani_2009")
