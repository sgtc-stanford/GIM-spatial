#!/bin/bash
# SG 6/18/2020: Adapted from: https://blog.sleeplessbeastie.eu/2019/11/11/how-to-parse-ini-configuration-file-using-bash/

# Get INI section
ReadINISections(){
  local filename="$1"
  awk '{ if ($1 ~ /^\[/) section=tolower(gensub(/\[(.+)\]/,"\\1",1,$1)); config[section]=1 } END {
         for (key in config) { print key} }' ${filename}
}

# Get/Set all INI sections
GetINISections(){
  local filename="$1"

  sections="$(ReadINISections $filename)"
  for section in $sections; do
    array_name="config_${section}"
    declare -g -A ${array_name}
  done
  eval $(awk -F= '{ 
                    if ($1 ~ /^\[/) 
                      section=tolower(gensub(/\[(.+)\]/,"\\1",1,$1)) 
                    else if ($1 !~ /^$/ && $1 !~ /^;/) {
                      gsub(/^[ \t]+|[ \t]+$/, "", $1); 
                      gsub(/[\[\]]/, "", $1);
                      gsub(/^[ \t]+|[ \t]+$/, "", $2); 
                      if (config[section][$1] == "")  
                        config[section][$1]=$2
                      else
                        config[section][$1]=config[section][$1]" "$2} 
                    } 
                    END {
                      for (section in config)    
                        for (key in config[section])  
                          print "config_"section"[\""key"\"]=\""config[section][key]"\";"
                    }' ${filename}
        )


}

# Parsing of input config file
if [ "$#" -eq "1" ] && [ -f "$1" ]; then
  filename="$1"
  GetINISections "$filename"

  for section in $(ReadINISections "$filename"); do
    echo "[${section}]"
    for key in $(eval echo $\{'!'config_${section}[@]\}); do
            echo -e "  ${key} = $(eval echo $\{config_${section}[$key]\})"
    done
  done
else
  echo "missing INI file"
fi
