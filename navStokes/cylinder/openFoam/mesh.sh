#!/usr/bin/bash

cd "${0%/*}" || exit 1

meshFile="cylinder.msh"

gmshToFoam $meshFile &&
  foamDictionary -set "entry0/empty/type=empty, entry0/wall/type=wall, entry0/cylinder/type=wall" \
    constant/polyMesh/boundary
