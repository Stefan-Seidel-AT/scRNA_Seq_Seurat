#! /bin/bash


echo "start Seurat :-) "
echo "$PWD"

### To start Seurat and mount current directory
#echo "docker run -it -v "$PWD":/work -p 8080:8080 satijalab/seurat:latest"
#docker run -it -v $PWD:/work -p 8080:8080 satijalab/seurat:latest
echo " to start R_docker insert: localhost:8787 to the browser "

#echo "docker run -v "$PWD":/work -p 8787:8787 -e DISABLE_AUTH=true 0f56500683f2"
echo "docker run -v "$PWD":/home/rstudio -p 8787:8787 -e DISABLE_AUTH=true 0f56500683f2"
#docker run -v $PWD:/work -p 8787:8787 -e DISABLE_AUTH=true 0f56500683f2
docker run -v $PWD:/home/rstudio -p 8787:8787 -e DISABLE_AUTH=true 0f56500683f2
