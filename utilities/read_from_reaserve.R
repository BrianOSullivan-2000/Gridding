
################################################################################

## Run commands in BASH to load NWP data from reaserve ##

################################################################################

read_from_reaserve <- function(var = "instant_uv_sfc_10", lead_time = "24",
                               year = "2024", month = "11", 
                               day = "01", hour = "00",
                               username = Sys.getenv("REASERVE_USERNAME"),
                               password = Sys.getenv("REASERVE_PASSWORD"), 
                               no_pwd = F, remote_path = NA){
  
  # var                   -variable of interest. I'm not fully sure about the
  #                       naming convention, but I found this link:
  # opendatadocs.dmi.govcloud.dk/Data/Forecast_Data_Weather_Model_HARMONIE
  #                       Note that I have not confirmed these are correct
  
  # year/month/day/hour   -specified time for model output
  # lead_time             -lead time for the forecast (i.e. how far into future)
  # username              -username for reaserve account
  # password              -password for reaserve account
  
  # no_pwd
  # I have my password saved as a local variable on my laptop (in .Renviron)
  # so I don't have to enter it every time. If you set no_pwd to zero, 
  # you'll get prompted to enter a password every time you call this function
  
  # remote_path           -if the templating with default arguments doesn't get
  #                       you the filepath you want, you can do it manually
  
  require(ssh); require(jsonlite)
  
  ## ssh into server
  if(no_pwd){
    session <- ssh_connect(paste0(username, "@reaserve"))
  }
  else{
    session <- ssh_connect(paste0(username, "@reaserve"), passwd = password)  
  }
  
  # Path of file you want - either using standard template or remote path
  if(is.na(remote_path)){
    remote_path <- paste0(
      "../../nfs/archive/prod/archive/Harmonie/UWCW_DINIeps/",
      year, "/", month, "/", day, "/", hour, "/mbr000/",
      "fc", year, month, day, hour, "00_", lead_time, "_", var, 
      "_geo_39.639_334.553_1909x1609x2000m.grib2")
  }
  
  res <- ssh_exec_internal(session, paste0("grib_ls -j ", remote_path))
  layers <- fromJSON(rawToChar(res$stdout))[[1]]$shortName
  
  for(n in 1:length(layers)){
    
    # Download from server
    cmd <- paste0(
      "grib_get_data -w shortName=", layers[n], " ", shQuote(remote_path),
      " | awk '{print $1\",\"$2\",\"$3}' ",
      "> temp_data", n, ".csv"
    )
    
    ssh_exec_wait(session, cmd)
    scp_download(session, paste0("temp_data", n, ".csv"), to = ".", verbose = F)
    
    if(n == 1){
      data = read.csv("temp_data1.csv")
      data[layers[1]] <- data$Value
      data$Value <- NULL
    }
    else{
      temp_data = read.csv(paste0("temp_data", n, ".csv"))
      data[layers[n]] = temp_data$Value
    }

    
    # Clean up
    ssh_exec_wait(session, paste0("rm -r temp_data", n, ".csv"))
    file.remove(paste0("temp_data", n, ".csv"))
  }

  ssh_disconnect(session)

  return(data)
}
