#!/bin/bash

hash svn 2>&- || { echo "0" ; exit 0 ; }
svn info | grep "vision" | head -n 1 | cut -d ":" -f 2