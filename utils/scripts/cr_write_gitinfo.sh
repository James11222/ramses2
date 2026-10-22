#!/bin/bash
#######################################
# cr_write_gitinfo.sh [output file] [patch dir] [build command]
# creates write_gitinfo.f90 with the git
# and build information (CMake build)
#######################################
GITBRANCH=$(git rev-parse --abbrev-ref HEAD)
GITHASH=$(git log --pretty=format:'%H' -n 1)
GITREPO=$(git config --get remote.origin.url)
BUILDDATE=$(date +"%D-%T")
PATCHDIR="$2"
BUILDCOMMAND=$(echo "$3" | sed 's/"/""/g')
cat <<EOF >$1
subroutine write_gitinfo
  use amr_parameters, ONLY:builddate,buildcommand,patchdir,gitrepo,gitbranch,githash

  builddate = "$BUILDDATE"
  buildcommand = "$BUILDCOMMAND"
  patchdir = "$PATCHDIR"
  gitrepo = "$GITREPO"
  gitbranch = "$GITBRANCH"
  githash = "$GITHASH"

  write(*,*)' '
  write(*,'(" compile date = ",A)')TRIM(builddate)
  write(*,'(" compile command = ",A)')TRIM(buildcommand)
  !  write(*,'(" patch dir    = ",A)')TRIM(patchdir)
  write(*,'(" remote repo = ",A)')TRIM(gitrepo)
  write(*,'(" local branch = ",A)')TRIM(gitbranch)
  write(*,'(" last commit = ",A)')TRIM(githash)
  write(*,*)' '

end subroutine write_gitinfo
EOF
