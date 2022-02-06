#!/bin/bash 


read -p 'Enter commit message: ' MESSAGE 

git add .

git commit -m "$MESSAGE"

git push


