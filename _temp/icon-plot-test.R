rm(list=ls())

library(grImport)                      # import eps files
setwd("./save/")
fragIcons <- PostScriptTrace("Frag-Black-icons.eps")
fragIcons <- readPicture("Frag-Black-icons.eps.xml")
setwd("..")
FragIcon1 <- fragIcons[49:50]          # 1. Continuous 
FragIcon2 <- fragIcons[29:48]          # 2. Corridors
FragIcon3 <- fragIcons[ 9:28]          # 3. Pseudo-corridors
FragIcon4 <- fragIcons[ 1:8 ]          # 4. Isolated

## str(fragIcons)
grid.picture(fragIcons)
grid.picture(fragIcons[1])             #  (first path)
grid.picture(fragIcons[1:8])           # 4. Isolated
grid.picture(fragIcons[9:28])          # 3. Pseudo-corridors
grid.picture(fragIcons[29:48])         # 2. Corridors
grid.picture(fragIcons[49:50])         # 1. Continuous

grid.picture(FragIcon1)
grid.picture(FragIcon2)
grid.picture(FragIcon3)
grid.picture(FragIcon4)

vp1 <- viewport(width = 0.1, height = 0.1, x = 1/6, y = -0.2)


