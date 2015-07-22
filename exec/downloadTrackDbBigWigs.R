argv <- "~/projects/blueprint/labels/H3K27ac_TDH_MonoMacroMyeloid/trackDb.txt"

trackDb.txt <- argv[1]

trackDb.lines <- readLines(trackDb.txt)
