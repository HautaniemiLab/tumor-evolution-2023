#!/bin/bash

#Short script for printing date, git version and full command to a log
#file at desired location

#Options
    #Defaults
    LOG_FILE="./$(date +"%Y%m%d.%H-%M-%S").log"
    COMMAND=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            --log-file )    shift
                            LOG_FILE=$1 ;;
            --command )     shift
                            COMMAND="$COMMAND $1" ;;
            * )             COMMAND="$COMMAND $1"
        esac
        shift
    done

#Make directory
mkdir -p $(dirname $LOG_FILE)

#Date and Git version
echo "#$(date)" > $LOG_FILE
echo "#$(git log -1 --oneline)" >> $LOG_FILE

#Uncommitted changes
git status -s | grep -v "^??" | sed 's/^/#/' >> $LOG_FILE

#Command itself with excess leading and trailing white spaces removed
echo "$COMMAND" | sed 's/^\s\+//; s/\s\+$//' >> $LOG_FILE

