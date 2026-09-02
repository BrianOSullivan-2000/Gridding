
################################################################################

## SQL Query function to load data in from sol ##

################################################################################

read_from_sol <- function(sql_query_string, 
                          username = Sys.getenv("SOL_USERNAME"),
                          password = Sys.getenv("SOL_PASSWORD"),
                          sol_direct = F){
  
  # sql_query_string      -string of SQL query for the data
  # username              -username for sol account
  # password              -password for sol account
  # sol_direct            -if already running on sol code is simpler
  
  # Note I'm using my own username/password saved as environment variables on
  # my laptop (.Renviron). For this function to work, you also need to be able 
  # to ssh into sol without a password, otherwise the code hangs and doesn't work.
  # This can be done by setting up a private key on your local laptop and
  # putting a matching public key on to sol.
  
  # Set sol_direct = TRUE if you are running this code directly on sol. Then
  # everything is much simpler using odbcConnect and sqlQuery directly
  
  require(RODBC); require(readr)
  
  if(!sol_direct){
    r_command <- paste0(
      "library(RODBC); ",
      "ch <- odbcConnect(\\\"climat\\\", \\\"", username,"\\\", \\\"", 
      password,"\\\"); ",
      "df <- sqlQuery(ch, \\\"", sql_query_string, "\\\"); ",
      "odbcClose(ch); ",
      "write.csv(df, file = '', row.names = TRUE, quote = TRUE)"
    )
    sol_command <- paste0(
      "ssh ", username, "@sol3 \"Rscript -e '", r_command, "'\""
    )
    
    csv_text <- system(sol_command, intern = TRUE)
    df <- read_csv(paste(csv_text, collapse = "\n"))
  }
  else{
    ch <- odbcConnect("climat", username, password)
    df <- sqlQuery(ch, sql_query_string)
    odbcClose(ch)
  }
  return(df)
}